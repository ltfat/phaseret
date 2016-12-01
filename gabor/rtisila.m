function [c,f,relres,iter]=rtisila(s,g,a,M,varargin)
%RTISILA Real-Time Iterative Spectrogram Inversion with Look Ahead
%   Usage: c = rtisila(s,g,a,M)
%          c = rtisila(s,g,a,M,maxit)
%          [c,f,relres,iter] = rtisila(...)
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
%   `rtisila(s,g,a,M)` attempts to find Gabor coefficients *c* such
%   that::
%
%     s = abs(c);
%
%   using the Real-Time Iterative Spectrogram Inversion with Look Ahead.
%
%   `[c,f,relres,iter]=rtisila(...)` additionally returns the final
%   residual `relres`, the number of per-frame iterations done `iter` and the
%   coefficients *c* with the reconstructed phase. The relationship between
%   *f* and *c* is::
%
%     f = idgtreal(c,gd,a,M)
%
%   where *gd* is the canonical dual window obtained by |gabdual|.
%
%   `rtisila` takes the following addition arguments:
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
%   See also:  lertisila, idgtreal, dgtreal
%
%   References: zhbewy07
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
[flags,kv,Ls]=ltfatarghelper({'maxit'},definput,varargin);
Ls = kv.Ls;

complainif_notposint(kv.maxit,'maxit',mfilename);

% Default number of lookahead frames.
if isempty(kv.lookahead)
    kv.lookahead = min([ceil(M/a)-1,N-1]);
end

if kv.lookahead < 0 || kv.lookahead > N - 1
    error('%s: lookahead must be in range [0-%d]',upper(mfilename),N-1);
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

[gnum,gdnum,specg1,specg2] = comp_rtisilawins(gnum,gd,a,M);

% Buffer initialization
cframes = zeros(M ,lookback + kv.lookahead+1);
sframes = zeros(M2,kv.lookahead+1);

c = zeros(M2,N);

% Preread modulus
sframes(:,2:end) = abss(:,1:kv.lookahead);

% n -th frame is the submit frame
for n=1:N
    % Index of the submit and the look-ahead frames in the spectrogram
    nextnewframeidx = mod(n - 1 + kv.lookahead,N) + 1;

    % Shift cols in the buffer
    cframes(:,1:end-1)   = cframes(:,2:end);
    % New empty column
    cframes(:,end) = 0;

    sframes(:,1:end-1)   = sframes(:,2:end);
    sframes(:,end) = abss(:,nextnewframeidx);

    % Update the lookahead frames and the submit frame 
    [cframes, c(:,n)] = ...
        comp_rtisilaupdate(cframes,gnum,specg1,specg2,gdnum,a,M,sframes,kv.lookahead,kv.maxit);

    % Get the submit frame
    % recframes(:,n) = cframes(:,lookback+1);
end

iter = kv.maxit*kv.lookahead;

if ~flags.do_timeinv
    c = phaseunlockreal(c,a,M);
end

f = idgtreal(c,gd,a,M,flags.phase);
% Cur or extend and reformat f
f = comp_sigreshape_post(f,Ls,0,[0; W]);

if nargout>2
    norm_s = norm(abss,'fro');
    relres = norm(dgtreal(f,g,a,M,flags.phase)-abss,'fro')/norm_s;
end








