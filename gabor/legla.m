function [c,f,relres,iter]=legla(s,g,a,M,varargin)
%LEGLA Le Roux's Griffin-Lim Algorithm for real signals
%   Usage: c = legla(s,g,a,M)
%          c = legla(s,g,a,M,maxit)
%          [c,f,relres,iter] = legla(...)
%
%   Input parameters:
%         s       : Array of initial coefficients.
%         g       : Analysis Gabor window
%         a       : Hop factor
%         M       : Number of channels
%         Ls      : length of signal.
%   Output parameters:
%         c       : Coefficients with the reconstructed phase
%         f       : Signal.
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%
%   `legla(s,g,a,M)` attempts to find coefficients *c* from their abs. 
%   value::
%
%     s = abs(dgtreal(f,g,a,M));
%
%   using Le Rouxs modifications of the Griffin-Lim algorithm.
%
%   `[c,f,relres,iter]=legla(...)` additionally returns an array
%   of residuals `relres`, the number of iterations done `iter` and the
%   reconstructed signal *f*. The relationship between *f* and *c* is::
%
%     f = idgtreal(c,gd,a,M)
%
%   where *gd* is the canonical dual window obtained by |gabdual|.
%
%   `legla` takes the following additional parameters:
%
%   Initial phase guess:
%
%     'input'      Choose the starting phase as the phase of the input
%                  *s*. This is the default
%
%     'zero'       Choose a starting phase of zero.
%
%     'rand'       Choose a random starting phase.
%
%   Enforcing prior information:
%
%     'coefmod',coefmod   Anonymous function in a form coefmod = @(c) ...;
%                         altering coefficients in each iteration after
%                         the phase update has been done.
%                         This is usefull when e.g. phase of some of
%                         the coefficients is known.
%
%   Projection kernel
%   -----------------
%   
%   The algorithm employs a twisted convolution of coefficients with the
%   truncated projection kernel. The full-size kernel is obtained as::
%
%     kern = dgt(gd,g,a,M)
%
%   where *gd* is canonical dual window obtained by |gabdual|. The
%   following key-value pairs control the final kernel size:
%
%     'relthr',relthr    The kernel is truncated such that it contains 
%                        coefficients with abs. values greater or equal 
%                        to *relthr* times the biggest coefficient. 
%                        The default value is 1e-3.
%
%     'kernsize',[height,width]  Define kernel size directly. When used,
%                                *relthr* is ignored.
%
%   Additinally, the phase update strategy is controlled by the following
%   flags:
%
%   Variant of the algorithm:
%
%     'trunc'      The truncated projection kernel is used directly.
%                  This is the default.
%
%     'modtrunc'   Modified phase update is done by setting the central
%                  sample of the projection kernel to zero.
%
%   Phase update scheme:
%
%     'stepwise'   Phase is updated only after the whole projection has
%                  been made.
%                  This is the default.
%
%     'onthefly'   The phase is updated for each coefficient immediatelly.
%
%   Algorithm acceleration:
%
%     'legla'      The original Giffin-Lim iteration scheme.
%                  This is the default.
%
%     'flegla'     A fast Griffin-Lim iteration scheme from Perraudin et. al.
%
%     'alpha',a    Parameter of the Fast Griffin-Lim algorithm. It is
%                  ignored if not used together with 'flegla' flag.
%                  The default value is 0.99.
%
%   Other:
%
%     'maxit',n    Do at most n iterations.
%
%     'print'      Display the progress.
%
%     'quiet'      Don't print anything, this is the default.
%
%     'printstep',p  If 'print' is specified, then print every p'th
%                    iteration. Default value is p=10;
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!
%   This function requires LTFAT 2.1.2 and above.
%   !!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!
%
%   References: leroux08 leroux10 pebaso13
%

%   AUTHOR: Zdenek Prusa

complainif_notenoughargs(nargin,4,mfilename);

if ~isnumeric(s)
    error('%s: s s must be numeric.',upper(mfilename));
end

definput.keyvals.Ls=[];
definput.keyvals.maxit=100;
definput.flags.startphase={'input','zero','rand'};
definput.flags.algvariant={'trunc','modtrunc'};
definput.flags.updatescheme={'stepwise','onthefly'};
definput.flags.method={'legla','flegla'};
definput.keyvals.alpha=0.99;
definput.flags.print={'quiet','print'};
definput.flags.phase={'timeinv','freqinv'};
definput.keyvals.relthr = 1e-3;
definput.keyvals.kernsize = [];
definput.keyvals.printstep=10;
definput.keyvals.coefmod = [];
[flags,kv]=ltfatarghelper({'maxit'},definput,varargin);
Ls = kv.Ls;

if ~isempty(kv.coefmod) && isa(kv.coefmod,'function_handle')
    error('%s: coefmod must be anonymous function.',upper(mfilename))
end

if kv.relthr>1 && kv.relthr<0
    error('%s: relthr must be in range [0,1].',upper(mfilename));
end

if flags.do_input
    % Start with the phase given by the input.
    c=s;
    if flags.do_timeinv
       c = phaseunlockreal(c,a,M);
    end
end;

if flags.do_zero
    % Start with a phase of zero.
    c=abs(s);
end;

if flags.do_rand
    c=abs(s).*exp(2*pi*1i*rand(size(s)));
end;

s = abs(s);

% For normalization purposes
norm_s=norm(s,'fro');

relres=zeros(kv.maxit,1);

[~,N,W] = size(s);
L = N*a;

% Projection kernel, this has the side effect that it checks inputs
projfnc = @(c) dgtreal(idgtreal(c,{'dual',g},a,M),g,a,M);
ctmp = zeros(floor(M/2)+1,N); ctmp(1) = 1;
kern = projfnc(ctmp);
clear ctmp;

% Replace projection with faster version
gnum = gabwin(g,a,M,L);
gd = gabdual(g,a,M,L);
projfnc = @(c) comp_sepdgtreal(comp_isepdgtreal(c,gd,L,a,M,0),gnum,a,M,0);

if isempty(kv.kernsize)
    kv.kernsize = findsmallkernelsize(kern,kv.relthr);
else
    if ~isnumeric(kv.kernsize) || numel(kv.kernsize)~=2
        error('%s: Kernel size must be a 2 element vector.',upper(mflename));
    end
%     if any(mod(kv.kernsize,2)==0)
%         error('%s: Kernel size must be odd.',upper(mfilename));
%     end
    if any(kv.kernsize<=0) || kv.kernsize(1)>M || kv.kernsize(2)>N
        error('%s: Invalid kernel size.',upper(mfilename));
    end
end

% Shrink kernel
kernsmall = postpad(kern,floor(kv.kernsize(1)/2)+1);
kernsmall = middlepad(kernsmall,kv.kernsize(2),2);

if flags.do_modtrunc
    kernsmall(1,1) = 0;
end

% Do explicit coefmod
if ~isempty(kv.coefmod)
    c = kv.coefmod(c);
end

if flags.do_legla
    for iter=1:kv.maxit

        if nargout>1
            % We do the projection explicitly as comp_leglaupdatereal
            % returns something different. This effectivelly doubles the
            % execution time.
            cproj = projfnc(c);
            relres(iter)=norm(abs(cproj)-s,'fro')/norm_s;
        end

        % Do the leGLA phase update
        c = comp_leglaupdate(c,kernsmall,s,a,M,flags.do_onthefly);
        
        % Apply coefficient restriction
        if ~isempty(kv.coefmod)
            c = kv.coefmod(c);
        end

        if flags.do_print
            if mod(iter,kv.printstep)==0
                fprintf('LEGLA: Iteration %i, residual = %f.\n',iter,relres(iter));
            end
        end
    end
elseif flags.do_flegla
    told=c;
    for iter=1:kv.maxit
        
        if nargout>1
            cproj = projfnc(c);
            relres(iter)=norm(abs(cproj)-s,'fro')/norm_s;
        end

        % Do the leGLA phase update
        tnew = comp_leglaupdate(c,kernsmall,s,a,M,flags.do_onthefly);
        
        % Apply coefficient restriction
        if ~isempty(kv.coefmod)
            tnew = kv.coefmod(tnew);
        end
        
        % The acceleration step
        c = tnew + kv.alpha*(tnew-told);
        
        % Keep for next iteration
        told = tnew;

        if flags.do_print
            if mod(iter,kv.printstep)==0
                fprintf('LEGLA: Iteration %i, residual = %f.\n',iter,relres(iter));
            end
        end
    end
end

f = idgtreal(c,gd,a,M,Ls);
f = comp_sigreshape_post(f,Ls,0,[0; W]);

if flags.do_timeinv
    c = phaselockreal(c,a,M);
end


function ksize=findsmallkernelsize(kern,relthr)
[M2,N] = size(kern);
thr = relthr*max(abs(kern(:)));

lastrow = 1;
for n=1:N
    newlastrow = find(abs(kern(:,n))>=thr,1,'last');
    if newlastrow > lastrow
        lastrow = newlastrow;
    end
end

lastcol1 = 1;
% Searching from the other end is not neccesary since the kenel
% is always symmetric in the horizontal direction too
% lastcol2 = 1;
for m=1:M2
    newlastcol = find(abs(kern(m,1:ceil(end/2)))>=thr,1,'last');
    if newlastcol > lastcol1
        lastcol1 = newlastcol;
    end
    
%     newlastcol = find(abs(kern(m,end:-1:floor(end/2)))>=thr,1,'last');
%     if newlastcol > lastcol2
%         lastcol2 = newlastcol;
%     end
end

ksize = [2*lastrow-1, 2*lastcol1 - 1];


