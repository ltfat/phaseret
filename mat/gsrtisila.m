function [c,f,relres,iter]=gsrtisila(s,g,a,M,varargin)
%GSRTISILA Gnann and Spiertz’s Real-Time Iterative Spectrogram Inversion
%   Usage: c = gsrtisila(s,g,a,M)
%          c = gsrtisila(s,g,a,M,Ls)
%          [c,f,relres,iter] = gsrtisila(...)
%
%   Input parameters:
%         s       : Modulus of coefficients.
%         g       : Analysis Gabor window
%         a       : Hop factor
%         M       : Number of channels
%         Ls      : length of signal.
%   Output parameters:
%         c       : Coefficients with the reconstructed phase.
%         f       : Reconstructed signal.
%         relres  : Final residual error.
%         iter    : Number of per-frame iterations done.
%
%   `gsrtisila(s,g,a,M)` attempts to find Gabor coefficients *c* of some
%   signal *f* such that::
%
%     c = dgtreal(f,g,a,M);
%     s = abs(c);
%
%   using the Gnann and Spiertz’s Real-Time Iterative Spectrogram
%   Inversion with Look Ahead.
%
%   `[c,f,relres,iter]=gsrtisila(...)` additionally returns the time domain
%   signal *f* and the residual error `relres`.
%   The relationship between *f* and *c* is the following::
%
%     f = idgtreal(c,gd,a,M)
%
%   where *gd* is the canonical dual window obtained by |gabdual|.
%
%   `gsrtisila` takes the following addition arguments:
%
%   Algorithm parameters:
%
%     'lookahead',lookahead  Number of lookahead frames. The default value
%                            is `ceil(M/a)-1`.
%
%     'maxit',maxit          Number of RTISILA iterations. The default
%                            value is 5. The total number
%                            of per-frame iteratins is `(lookahead+1)*maxit`.
%
%
%   See also:  idgtreal, dgtreal
%
%   References: gnsp08 gnsp10
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
definput.keyvals.lookahead = [];
definput.flags.phase={'freqinv','timeinv'};
definput.flags.lastbufinit={'zeros','unwrap','input','spsi','rtpghi'};
definput.keyvals.rtpghi = {};
definput.keyvals.unwrappar=0.3;
definput.flags.frameorder={'reverse','energy'};
[flags,kv,Ls]=ltfatarghelper({'Ls','maxit'},definput,varargin);

complainif_notposint(kv.maxit,'maxit',mfilename);

% Default number of lookahead frames.
if isempty(kv.lookahead)
    kv.lookahead = min([ceil(M/a)-1,N-1]);
end

if kv.lookahead < 0 || kv.lookahead > N - 1
    error('%s: lookahead must be in range [0-%d]',upper(mfilename),N-1);
end

if flags.do_rtpghi
    if isempty(kv.rtpghi)
        error('%s: RTPGHI parameters cell is empty.',upper(mfilename));
    end
    
    if ~iscell(kv.rtpghi) && isscalar(kv.rtpghi)
        kv.rtpghi = { kv.rtpghi };
    end
    
    definput2.keyvals.gamma = [];
    definput2.keyvals.tol=[1e-6];
    definput2.flags.variant={'normal','causal'};
    [flags2,~,gamma,tol]=ltfatarghelper({'gamma','tol'},definput2,kv.rtpghi);   
    if isempty(gamma)
        error('%s: RTPGHI gamma parameter is missing.',upper(mfilename));
    end
end


% Analysis window (as array)
gnum = gabwin(g,a,M);

if numel(gnum) > M
    error(['%s: The algorithm does not work for non-painless Gabor ',...
    'system, that is when numel(g)>M.'],upper(mfilename))
end

abss = abs(s);


lookback = max([ceil(M/a) - 1, kv.lookahead]);

% Synthesis window
gd = gabdual(g,a,M,L);
% ... as array
gdnum = gabwin(gd,a,M);
% Create and overlay several scaled dual windows
wins = repmat(gnum.*gdnum*M,1,lookback+1+kv.lookahead);
winsum = overlayframes(wins,a,M);
rellim = 1e-3;
idx = abs(winsum)<rellim;
winsum(idx & winsum > 0) = rellim;
winsum(idx & winsum == 0) = 1;
winsum(idx & winsum < 0) = -rellim;

gnums = zeros(numel(gnum),kv.lookahead+1);
for la = 1:kv.lookahead + 1
    idx = la + lookback - 1;
    gnums(:,la) = fftshift(gnum)./winsum(a*idx+1:a*idx + M);
end

gdnum = fftshift(gdnum);

% Buffer initialization
cframes = zeros(M ,lookback + kv.lookahead+1);
sframes = zeros(M2,kv.lookahead+1);
coefbuf = zeros(M2 ,lookback + kv.lookahead+1);

if flags.do_rtpghi
    [tgrad, fgrad] = comp_pghiphasegrad( abss, gamma, a, M, 1, flags2.do_causal);
    tmpmask = zeros(M2,2); tmpmask(:,1) = 1;
end


c = zeros(M2,N);
% Preread modulus
sframes(:,2:end) = abss(:,1:kv.lookahead);

% n -th frame is the submit frame
for n=1:N
    % Index of the submit and the look-ahead frames in the spectrogram
    nextnewframeidx = mod(n - 1 + kv.lookahead,N) + 1;
    nextframeidx = mod(n - 1 - 1 + kv.lookahead,N) + 1;
    idx = mod(n - 1 - 1 + kv.lookahead:n - 1 + kv.lookahead,N) + 1;

    sframes(:,1:end-1)   = sframes(:,2:end);
    sframes(:,end) = abss(:,nextnewframeidx);

    % Shift cols in the buffer
    cframes(:,1:end-1)   = cframes(:,2:end);
    cframes(:,end) = 0;

    % Coefbuf
    coefbuf(:,1:end-1)   = coefbuf(:,2:end);
    coefbuf(:,end) = 0;

    if flags.do_spsi
        oldphase = angle(coefbuf(:,end-1));
        coefbuf(:,end) = comp_spsi(sframes(:,end),a,M,oldphase);
        cframes(:,end) = gdnum.*fftshift(comp_ifftreal(coefbuf(:,end),M))*M;
    elseif flags.do_unwrap
        nom = sframes(:,end).*coefbuf(:,end-1).^2.*sframes(:,end-2);
        denom = sframes(:,end-1).^2.*coefbuf(:,end-2);
        coefbuf(:,end) = nom./ denom;
        coefbuf(abs(denom)<1e-8,end) = 0;
        cframes(:,end) = 0.0*fftshift(comp_ifftreal(coefbuf(:,end),M))*M;
    elseif flags.do_input
        cframes(:,end) = gdnum.*fftshift(comp_ifftreal(s(:,nextnewframeidx),M))*M;
    elseif flags.do_rtpghi
        oldphase = angle(coefbuf(:,end-1:end));
        newphasetmp = comp_constructphasereal(abss(:,[nextframeidx,nextnewframeidx]),tgrad(:,idx),fgrad(:,idx),a,M,tol,2,tmpmask,oldphase);
        coefbuf(:,end) = sframes(:,end).*exp(1i*newphasetmp(:,2));
        cframes(:,end) = gdnum.*fftshift(comp_ifftreal(coefbuf(:,end),M))*M;
    end

    % Update the lookahead frames and the submit frame
    [cframes, coefbuf, c(:,n)] = ...
    comp_gsrtisilaupdate(cframes,coefbuf,gnums,gdnum,a,M,sframes,kv.lookahead,kv.maxit,flags.do_energy);
end

iter = kv.maxit*kv.lookahead;

% Alternative way of reconstruction from recframes
%
% f = zeros(L,1);
% idxrange = [0:floor(M/2),-ceil(M/2)+1:-1];
% for n=0:N-1
%     idx = mod(n*a + idxrange,L) + 1;
%     f(idx) = f(idx) + recframes(:,n+1);
% end
% c = dgtreal(f,g,a,M,flags.phase);

if ~flags.do_timeinv
    c = phaseunlockreal(c,a,M);
end

f = idgtreal(c,gd,a,M,flags.phase);

if nargout > 2
    norm_s = norm(abss,'fro');
    relres = norm(dgtreal(f,g,a,M,flags.phase)-abss,'fro')/norm_s;
end

% Cur or extend and reformat f
f = comp_sigreshape_post(f,Ls,0,[0; W]);

function partrec = overlayframes(cframes,a,M)

N = size(cframes,2);
bufLen = N*a - (a-1) + M-1;
partrec = zeros(bufLen,1);

startidx = ceil(M/2)-1;
idxrange = startidx + [0:floor(M/2),-ceil(M/2)+1:-1];
for n=0:N-1
    idx = n*a + idxrange + 1;
    partrec(idx) = partrec(idx) + cframes(:,n+1);
end






