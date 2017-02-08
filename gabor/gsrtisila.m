function [c,f,relres,iter]=gsrtisila(s,g,a,M,varargin)
%GSRTISILA Gnann and Spiertz’s Real-Time Iterative Spectrogram Inversion
%   Usage: c = gsrtisila(s,g,a,M)
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
%     c = dgtreal(f,g,a,M,'timeinv');
%     s = abs(c);
%
%   using the Gnann and Spiertz’s Real-Time Iterative Spectrogram
%   Inversion with Look Ahead.
%
%   `[c,f,relres,iter]=gsrtisila(...)` additionally returns the time domain
%   signal *f* and the residual error `relres`.
%   The relationship between *f* and *c* is the following::
%
%     f = idgtreal(c,gd,a,M,'timeinv')
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
%     'Ls',Ls                Crop the reconstructed signal *f* to length *Ls*.
%
%   Phase initialization:
%   ---------------------
%
%   Optionally, the following phase initialization strategies for the
%   newest lookahead frame can be employed:
%
%      'zeros'          Initialize with zeros. This is the default.
%
%      'input'          Use phase of *s*.
%
%      'unwrap'         Use phase vocoder style phase unwrapping.
%
%      'spsi'           Use the SPSI algorithm.
%
%      'rtpghi',gamma   Use the RTPGHI algorithm. The algorithm requires
%                       one additional look-ahead frame.
%
%      'rtpghi',{gamma,'causal'}  Use the causal version of RTPGHI.
%
%      'rtpghi',{gamma,tol,...}   Set tolerance of RTPGHI to `tol`. 
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
definput.flags.phase={'timeinv','freqinv'};
definput.flags.lastbufinit={'zeros','unwrap','input','spsi','rtpghi'};
definput.keyvals.rtpghi = {};
definput.keyvals.unwrappar=0.3;
definput.keyvals.mask=[];
definput.keyvals.usephase=[];
[flags,kv]=ltfatarghelper({'mask','usephase'},definput,varargin);
Ls = kv.Ls;

if ~isempty(kv.mask) || ~isempty(kv.usephase)
    error('%s: TODO: ''mask'' or ''usephase'' parameters are not supported yet.',...
          upper(mfilename));
end

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
    definput2.keyvals.tol=1e-6;
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
[gnums,gdnum] = comp_gsrtisilawins(gnum,gd,a,M,kv.lookahead);

% Buffer initialization
frames = zeros(M ,lookback + kv.lookahead+1);
sframes = zeros(M2,kv.lookahead+1);
cframes = zeros(M2 ,lookback + kv.lookahead+1);

if flags.do_rtpghi
    [tgrad, fgrad, logs] = comp_pghiphasegrad( abss, gamma, a, M, 1, flags2.do_causal);
end

c = zeros(M2,N);
% Preread modulus
sframes(:,2:end) = abss(:,1:kv.lookahead);

if flags.do_input
    cframes(:,2:end) = s(:,1:kv.lookahead);
end

if flags.do_unwrap
   omega = 2*pi*a*(0:M2-1)'/M;
end

% n -th frame is the submit frame
for n=1:N
    % Index of the submit and the look-ahead frames in the spectrogram
    nextnewframeidx = mod(n - 1 + kv.lookahead,N) + 1;
    idx = mod(n - 1 - 1 + kv.lookahead:n - 1 + kv.lookahead,N) + 1;

    [frames,cframes,sframes] = shiftcolsleft(frames,cframes,sframes);
    sframes(:,end) = abss(:,nextnewframeidx);

    if flags.do_spsi
        oldphase = angle(cframes(:,end-1));
        cframes(:,end) = comp_spsi(sframes(:,end),a,M,oldphase);
    elseif flags.do_unwrap
        phase2 = angle(cframes(:,end-2));
        phase1 = angle(cframes(:,end-1));
        phase0 = phase1 + omega + vocoderprincarg(phase1 - phase2 - omega);
        cframes(:,end) = kv.unwrappar*sframes(:,end).*exp(1i*phase0);
    elseif flags.do_input
        cframes(:,end) = s(:,nextnewframeidx);
    elseif flags.do_rtpghi
        oldphase = angle(cframes(:,end-1));
        newphase = comp_rtpghiupdate(logs(:,idx),tgrad(:,idx),fgrad(:,idx(2)),oldphase,tol,M);
        cframes(:,end) = sframes(:,end).*exp(1i*newphase);
    end

    % Update the lookahead frames and the submit frame
    [frames, cframes, c(:,n)] = ...
    comp_gsrtisilaupdate(frames,cframes,gnums,gdnum,a,M,sframes,kv.lookahead,kv.maxit);
end

iter = kv.maxit*kv.lookahead;

if ~flags.do_timeinv
    c = phaseunlockreal(c,a,M);
end

f = idgtreal(c,gd,a,M,Ls,flags.phase);

if nargout > 2
    norm_s = norm(abss,'fro');
    relres = norm(dgtreal(f,g,a,M,flags.phase)-abss,'fro')/norm_s;
end

function phase_out = vocoderprincarg(phase)
phase_out = phase - 2*pi*round(phase/(2*pi));

if any(abs(exp(1i*phase)-exp(1i*phase_out))>1e-10)
    error('Unwrapping failed');
end

function [frames,cframes,sframes] = shiftcolsleft(frames,cframes,sframes)

sframes(:,1:end-1)   = sframes(:,2:end);
sframes(:,end) = 0;

% Shift cols in the buffer
frames(:,1:end-1)   = frames(:,2:end);
frames(:,end) = 0;

% Coefbuf
cframes(:,1:end-1)   = cframes(:,2:end);
cframes(:,end) = 0;
