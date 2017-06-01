function [c,f,relres,iter]=gla(s,g,a,M,varargin)
%GLA Griffin-Lim Algorithm
%   Usage: c = gla(s,g,a,M)
%          c = gla(s,g,a,M,maxit)
%          c = gla(s,g,a,M,maxit,tol)
%          [c,f,relres,iter] = gla(...)
%
%   Input parameters:
%         s       : Initial coefficients.
%         g       : Analysis Gabor window
%         a       : Hop factor
%         M       : Number of channels
%         maxit   : Maximum number of iterations.
%         tol     : relative tolerance
%   Output parameters:
%         c       : Coefficients with the reconstructed phase
%         f       : Reconstructed signal.
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%
%   `gla(s,g,a,M)` attempts to find coefficients *c* from
%   their abs. value::
%
%     c = dgtreal(f,g,a,M,'timeinv');
%     s = abs(c);  
%
%   using the Griffin-Lim algorithm.
%
%   `[c,f,relres,iter]=gla(...)` additionally returns an array
%   of residuals `relres`, the number of iterations done `iter` and the
%   coefficients *c* with the reconstructed phase. The relationship between
%   *f* and *c* is::
%
%     f = idgtreal(c,gd,a,M,'timeinv')
%
%   where *gd* is the canonical dual window obtained by |gabdual|.
%
%   Initial phase guess
%   -------------------
%
%     'input'      Choose the starting phase as the phase of the input
%                  *s*. This is the default
%
%     'zero'       Choose a starting phase of zero.
%
%     'rand'       Choose a random starting phase.
%
%   Enforcing prior information
%   ---------------------------
%
%     'coefmod',coefmod   Anonymous function in a form coefmod = @(c) ...;
%                         altering coefficients in each iteration after
%                         the phase update has been done.
%                         This is usefull when e.g. phase of some of
%                         the coefficients is known.
%
%     'timemod',timemod   Anonymous function in a form timemod = @(f) ...;
%                         altering the time-domain signal in each iteration.
%                         This is usefull when e.g. the time support of the
%                         signal is known.
%                         Note that `numel(f)= size(s,2)*a`.
%
%   Algorithm acceleration
%   ----------------------
%
%     'gla'      The original Giffin-Lim iteration scheme.
%                This is the default.
%
%     'fgla'     A fast Griffin-Lim iteration scheme from Perraudin et. al..
%
%     'alpha',a    Parameter of the Fast Griffin-Lim algorithm. It is
%                  ignored if not used together with 'fgla' flag.
%
%   Additional parameters
%   ---------------------
%
%     'maxit',n    Do at most n iterations.
%
%     'tol',t      Stop if relative residual error is less than the
%                  specified tolerance.
%
%     'Ls',Ls      Crop the reconstructed signal *f* to length *Ls*.
%
%     'print'      Display the progress. This is disabled by default.
%
%     'printstep',p  If 'print' is specified, then print every p'th
%                    iteration. Default value is p=10;
%
%   See also: dgtreal idgtreal gabdual
%
%   References: griflim84 pebaso13
%

%   AUTHOR: Zdenek Prusa, Peter Soendergaard

definput.keyvals.Ls=[];
definput.keyvals.tol=1e-6;
definput.keyvals.maxit=100;
definput.flags.startphase={'input','zero','rand'};
definput.flags.method={'gla','fgla'};
definput.keyvals.alpha=0.99;
definput.flags.print={'quiet','print'};
definput.flags.phase={'timeinv','freqinv'};
definput.keyvals.kernsize = [];
definput.keyvals.printstep=10;
definput.keyvals.coefmod = [];
definput.keyvals.timemod = [];
[flags,kv]=ltfatarghelper({'maxit','tol'},definput,varargin);
Ls = kv.Ls;

if ~isempty(kv.coefmod) && isa(kv.coefmod,'function_handle')
    error('%s: coefmod must be anonymous function.',upper(mfilename))
end

if ~isempty(kv.timemod) && isa(kv.timemod,'function_handle')
    error('%s: timemod must be anonymous function.',upper(mfilename))
end

[M2,N,W] = size(s);
L = N*a;

M2true = floor(M/2) + 1;

if M2true ~= M2
    error('%s: Mismatch between *M* and the size of *s*.',thismfilename);
end

if flags.do_input
    % Start with the phase given by the input.
    c=s;
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

gnum = gabwin(g,a,M,L);
gd = gabdual(g,a,M,L);

fwdtra = @(f)  comp_sepdgtreal(f,gnum,a,M,flags.do_timeinv);
backtra = @(c) comp_isepdgtreal(c,gd,L,a,M,flags.do_timeinv);

% Do explicit coefmod
if ~isempty(kv.coefmod)
    c = kv.coefmod(c);
end

if flags.do_gla
    for iter=1:kv.maxit
        fiter = backtra(c);

        if ~isempty(kv.timemod)
            fiter = kv.timemod(fiter);
        end

        c = fwdtra(fiter);

        relres(iter) = norm(abs(c)-s,'fro')/norm_s;

        c = s.*exp(1i*angle(c));

        if ~isempty(kv.coefmod)
            c = kv.coefmod(c);
        end

        if relres(iter)<kv.tol
            relres=relres(1:iter);
            break;
        end;

        if flags.do_print
            if mod(iter,kv.printstep)==0
                fprintf('LEGLA: Iteration %i, residual = %f.\n',iter,relres(iter));
            end
        end
    end
elseif flags.do_fgla
    told=c;
    for iter=1:kv.maxit
        % Synthesis
        fiter = backtra(c);

        % Put restriction on f
        if ~isempty(kv.timemod)
            fiter = kv.timemod(fiter);
        end

        c = fwdtra(fiter);

        relres(iter)=norm(abs(c)-s,'fro')/norm_s;

        % Phase update
        tnew = s.*exp(1i*angle(c));

        % The acceleration step
        c = tnew + kv.alpha*(tnew-told);

        % Put restriction on c
        if ~isempty(kv.coefmod)
            c = kv.coefmod(c);
        end

        told = tnew;

        if relres(iter)<kv.tol
            relres=relres(1:iter);
            break;
        end;

        if flags.do_print
            if mod(iter,kv.printstep)==0
                fprintf('GLA: Iteration %i, residual = %f.\n',iter,relres(iter));
            end
        end

    end
end

f = idgtreal(c,gd,a,M,Ls,flags.phase);
f = comp_sigreshape_post(f,Ls,0,[0; W]);


