function [c,relres,iter,f]=lertisila(s,g,a,M,varargin)
%LERTISILA RTISI for real signals using Le Roux's truncated updates
%   Usage: c = lertisila(s,g,a,M)
%          c = lertisila(s,g,a,M,Ls)
%          [c,relres,iter,f] = lertisila(...)
%
%   Input parameters:
%         s       : Array of initial coefficients.
%         g       : Analysis Gabor window
%         a       : Hop factor
%         M       : Number of channels
%         Ls      : length of signal.
%   Output parameters:
%         f       : Reconstructed signal.
%         relres  : Final residual error.
%         iter    : Number of per-frame iterations done.
%         c       : Coefficients with the reconstructed phase.
%
%   `lertisila(s,g,a,M)` attempts to find Gabor coefficients *c* such
%   that::
%
%     s = abs(c);
%
%   using the Real-Time Iterative Spectrogram Inversion with Look Ahead
%   and Le Roux's truncated phase updates.
%   
%   `[c,relres,iter,f]=lertisila(...)` additionally returns an array
%   of residuals `relres`, the number of per-frame iterations done `iter` 
%   and the coefficients *c* with the reconstructed phase. The relationship
%   between *f* and *c* is::
%
%     f = idgtreal(c,gd,a,M)
%
%   where *gd* is the canonical dual window obtained by |gabdual|.
%
%   `lertisila` takes the following addition arguments:
%
%   Algorithm parameters:
%
%     'lookahead',la   Number of lookahead frames. 
%                      The default value is `ceil(M/a)-1`.
%
%     'kernsize',ks    The truncated dimensions [heigh,width] of the kernel.
%                      The default value is `2*lookahead+1` for both.
%
%     'asymwin'        Use asymetric window for the newest lookahead
%                      frame.
%                      This is the default.
%
%     'regwin'         Use regular window for the newest lookahead
%                      frame.
%
%   Initial phase guess:
%
%     'zhu'        The initial phase of the newly read look-ahead frame
%                  is computed from the overlapped previous frames.
%                  This is the default.
%
%     'input'      Choose the starting phase as the phase of the input
%                  *s*.
%
%     'zero'       Choose a starting phase of zero.
%
%     'rand'       Choose a random starting phase.
%
%     'unwrap'     Determines phase by unwrapping the phase from the
%                  previous frames.
%
%     'unwrappar',upar  Argument of the unwrapping phase initialization.
%                       It is ignored if the 'unwrap' parameter is not
%                       used. 
%                       The default value is 0.3.
%
%   Variant of the algorithm:
%
%     'trunc'      The projection kernel is used directly.
%                  This is the default.
%
%     'modtrunc'   Modified phase update is done by setting the central
%                  sample of the projection kernel to zero.
%
%   Phase update scheme:
%
%     'framewise'  Phase is updated for each frame.
%                  This is the default.
%
%     'onthefly'   The phase is updated for each coefficient immediatelly.
%
%   Parameters of modifications:
%
%     'energy'     Process the lookahead frames in the order of their
%                  energy.
%   
%   References: zhbewy07 gnsp10 leroux10
%

%   AUTHORS: Zdenek Prusa   
%

complainif_notposint(a,'a',mfilename);
complainif_notposint(M,'M',mfilename);

if ~isnumeric(s) || isempty(s)
    error('%s: s must be a numeric array of coefficients.',upper(mfilename));
end

if size(s,3)>1
   error('%s: s cannot be 3dimensional.',upper(mfilename)); 
end

[M2,N,W] = size(s);
L = N*a;

definput.keyvals.Ls=[];
definput.keyvals.maxit=5;
definput.keyvals.unwrappar=0.3;
definput.flags.startphase={'zhu','input','zero','rand','unwrap'};
definput.flags.frameorder={'plain','energy'};
definput.flags.asymwin={'asymwin','regwin'};
definput.keyvals.lookahead = [];
definput.keyvals.kernsize = [];
definput.flags.phase={'freqinv','timeinv'};
definput.flags.algvariant={'trunc','modtrunc'};
definput.flags.updatescheme={'framewise','onthefly'};
[flags,kv,Ls]=ltfatarghelper({'Ls','maxit'},definput,varargin);

if ~isnumeric(kv.maxit) || any(kv.maxit<=0) || any(rem(kv.maxit,1))
    error('%s: maxit must be a positive integer (or array).',upper(mfilename));
end

% Default number of lookahead frames.
% The kernel size in the horizontal direction is 2*kv.lookahead + 1
if isempty(kv.lookahead)
    kv.lookahead = ceil(M/a)-1;
else
    if kv.lookahead < 0 || kv.lookahead > N - 1
        error('%s: lookahead must be in range [0-%d]',upper(mfilename),N-1);
    end
end

% Default kernel size in the frequency direction
if isempty(kv.kernsize)
    kv.kernsize = [2*(kv.lookahead) + 1, 2*(kv.lookahead) + 1];
end

% Kernel dimensions
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
kernh = kv.kernsize(1);
kernw = kv.kernsize(2);
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if kernh < 1 || kernh > M2
   error(['%s: Kernel size in the frequency direction must be',...
          ' in range [1-%d]'],upper(mfilename),M2); 
end

if kernw > 2*(kv.lookahead) + 1
    error(['%s: Kernel size in the time dimension must not be bigger',...
           ' than 2*(lookahead) + 1'],upper(mfilename));
end

if any(mod(kv.kernsize,2)==0)
    error('%s: Kernel size must be odd.',upper(mfilename));
end

abss = abs(s);
norm_s = norm(abss,'fro');
if flags.do_input
  % Start with the phase given by the input.
  c=s;
end;

if flags.do_zero || flags.do_unwrap || flags.do_zhu
  % Start with a phase of zero.
  c=abss;
end;

if flags.do_rand
  % Start with a renadom phase
  c=abss.*exp(2*pi*1i*rand(size(s)));
end;

% Dual window
gd = gabdual(g,a,M,L);

% Prepare truncated projection kernels
projfncBase = @(c) dgt(idgt(c,gd,a),g,a,M);
projfncSpec = @(c) dgt(idgt(c,gd,a),gd,a,M);
projfncBaseReal = @(c) dgtreal(idgtreal(c,gd,a,M),g,a,M);
ctmp = zeros(M,N); ctmp(1) = 1;
kern = projfncBase(ctmp);

if flags.do_modtrunc
    kern(1,1) = 0;
end

% Just shrink the kernel to the size of look ahead
kernsmall = middlepad2(kern,[kernh, kernw]);

% Kernel for the first asymetric window
ctmpSpec1 = zeros(M,N);
ctmpSpec1(1,[2:kernw+1]) = 1;
kernSpec1 = projfncSpec(ctmpSpec1);

if flags.do_modtrunc
    kernSpec1(1,1) = 0;
end

% Just shrink the kernel to the size of look ahead
kernsmallSpec1 = middlepad2(kernSpec1,[kernh, kernw]);

% Kernel for the second asymetric window
ctmpSpec2 = zeros(M,N);
ctmpSpec2(1,[1:kernw+1]) = 1;
kernSpec2 = projfncSpec(ctmpSpec2);

if flags.do_modtrunc
    kernSpec2(1,1) = 0;
end

% Just shrink the kernel to the size of look ahead
kernsmallSpec2 = middlepad2(kernSpec2,[kernh,kernw]);

cfull = c;
c = zeros(size(cfull));

cbuf = zeros(M2,3*(kv.lookahead) + 1);

% Precompute modulated kernels
kNo = lcm(M,a)/a;
kernelsSmall = zeros([kernh,kernw,kNo]);
kernelsSpec1Small = kernelsSmall;
kernelsSpec2Small = kernelsSmall;
for k = 0:kNo-1
    kernelsSmall(:,:,k+1) = phasekernfi(kernsmall,k,a,M);
    kernelsSpec1Small(:,:,k+1) = phasekernfi(kernsmallSpec1,k,a,M);
    kernelsSpec2Small(:,:,k+1) = phasekernfi(kernsmallSpec2,k,a,M);
    
    kernelsSmall(:,:,k+1)=fftshift(involute2(conj(kernelsSmall(:,:,k+1))));
    kernelsSpec1Small(:,:,k+1)=fftshift(involute2(conj(kernelsSpec1Small(:,:,k+1))));
    kernelsSpec2Small(:,:,k+1)=fftshift(involute2(conj(kernelsSpec2Small(:,:,k+1))));
end

% Buffer initialization
cbuf(:,1:2*kv.lookahead+1) = cfull(:,mod( -1 + (-kv.lookahead:kv.lookahead), N)+1);

kernw2 = floor(kernw/2);

% n -th frame is the submit frame
for n=1:N
    % Shift cols in the buffer 
    cbuf(:,1:kernw-1) = cbuf(:,2:kernw);
    cbuf(:,kernw:end) = 0;
    
    % Index of the newest look-ahead frame
    nextnewframeidx = mod(n - 1 + kv.lookahead,N) + 1;
    
    % Pick new lookahead frame
    if flags.do_unwrap
        newframeidx1 = mod(nextnewframeidx-1 -1,N) +1;
        newframeidx2 = mod(nextnewframeidx-2 -1,N) +1;
        
        cbuf(:,kernw) = kv.unwrappar*abs(cfull(:,nextnewframeidx)) ...
                                    .*cfull(:,newframeidx1).^2 ...
                                    .*abs(cfull(:,newframeidx2))./ ...
                                    (abs(cfull(:,newframeidx1)).^2 ...
                                    .*cfull(:,newframeidx2) );
    elseif ~flags.do_zhu
        % Use the initial estimate
        cbuf(:,kernw) = cfull(:,nextnewframeidx);      
    elseif flags.do_zhu
        % Do nothing. The original algorithm assumes the new frame to only
        % consist of the sum of the overlapping preceeding frames.
        % This is equivalent to setting cbuf(:,kernw) to zeros.
        % It gets filled after the first iteration.
    end
    
    %% Explicit first iteration
    % 1) No energy sorting
    % 2) Special analysis window for the rightmost frame 
    %    There are two different "special" analysis frames
    %    depending on flags.do_zhu
    
    %% 1) The new looakahead frame
    if flags.do_asymwin
        if flags.do_zhu
           kernact = kernelsSpec1Small(:,:,mod(n-1 + kv.lookahead,kNo)+1);
        else
           kernact = kernelsSpec2Small(:,:,mod(n-1 + kv.lookahead,kNo)+1);
        end
    else
        kernact = kernelsSmall(:,:,mod(n-1 + kv.lookahead,kNo)+1);
    end

    indx = kv.lookahead+kv.lookahead+1;
    indxRange = indx + (-kernw2:kernw2);
    
    cbuf(:,indx) = comp_leglaupdatesinglecol(cbuf(:,indxRange),...
                         kernact,abs(cfull(:,nextnewframeidx)),M,flags.do_onthefly);

    
    %% 2) Other lookahead frames and the submit frame
    for nback = kv.lookahead-1:-1:0
       kernact = kernelsSmall(:,:,mod(n-1 + nback,kNo)+1);
       indx = kv.lookahead+nback+1;
       indxRange = indx + (-kernw2:kernw2);
       
       cbuf(:,indx) = comp_leglaupdatesinglecol(cbuf(:,indxRange),...
                         kernact,abs(cfull(:,mod(n-1+nback,N)+1)),M,flags.do_onthefly);
    end
    
    % Determine order of frames for the next iterations
    frameorder = kv.lookahead:-1:0;
    if flags.do_energy
        [~,frameorder] = sort(sum(abs(cbuf(:,kv.lookahead+1:kernw)).^2),'descend');
        frameorder = frameorder - 1;
    end
      
    %% 3) All the other iterations
    for iter=2:kv.maxit
        for nback = frameorder
            % Last frame gets (another) special analysis window
            if flags.do_asymwin && nback == kv.lookahead
                kernact = kernelsSpec2Small(:,:,mod(n-1 + kv.lookahead,kNo)+1);
            else
                % Pick the right kernel
                kernact = kernelsSmall(:,:,mod(n-1 + nback,kNo)+1);
            end
            
            indx = kv.lookahead+nback+1;
            indxRange = indx + (-kernw2:kernw2);
            
            cbuf(:,indx) = comp_leglaupdatesinglecol(cbuf(:,indxRange),...
                              kernact,abs(cfull(:,mod(n-1+nback,N)+1)),M,flags.do_onthefly);
        end
    end
    
    % Submit a frame
    c(:,n) = cbuf(:,kv.lookahead + 1);
end

relres = norm(abs(projfncBaseReal(c))-abss,'fro')/norm_s;
iter = kv.maxit*kv.lookahead;

f = idgtreal(c,gd,a,M,Ls);
f = comp_sigreshape_post(f,Ls,0,[0; W]);

if flags.do_timeinv
    % Convert to time invariant phase
    N = size(c,2);
    L = N*a;
    b = L/M;
    M2=floor(M/2)+1;
    N=size(chat,2);

    TimeInd = (0:(N-1))/N;
    FreqInd = (0:(M2-1))*b;

    phase = FreqInd'*TimeInd;
    phase = exp(2*1i*pi*phase);

    % Handle multisignals
    c = bsxfun(@times,c,phase);
end

% M/a periodic in n
function kernm = phasekernfi(kern,n,a,M)
kern = involute2(kern);

kernh = size(kern,1);
l = -2*pi*n*fftindex(kernh)*a/M;

kernm = bsxfun(@times,kern,exp(1i*l));

% M/a periodic in m
function kernm = phasekernti(kern,m,a,M)

[kernh, kernw] = size(kern);
l = 2*pi*m*fftindex(kernw)'*a/M;

kernm = bsxfun(@times,kern,exp(1i*l));


function f=involute2(f)
f = involute(f,1);
f = conj(involute(f,2));

function f=middlepad2(f,L,varargin)
f = middlepad(f,L(1),1,varargin{:});
f = middlepad(f,L(2),2,varargin{:});










