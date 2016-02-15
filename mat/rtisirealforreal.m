function [c,relres,iter,f]=rtisirealforreal(s,g,a,M,varargin)
%RTISIREAL Real-Time Iterative Spectrogram Inversion (RTISI) for real signals
%   Usage: f = rtisireal(s,g,a,M)
%          f = rtisireal(s,g,a,M,Ls)
%          [f,relres,iter,c] = rtisireal(...)
%
%   Input parameters:
%         s       : Modulus of coefficients.
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
%   `rtisireal(s,g,a,M)` attempts to find a signal *f* with has *s* as
%   the abs. value of the Gabor coefficients such as::
%
%     s = abs(dgtreal(f,g,a,M));
%
%   using the Real-Time Iterative Spectrogram Inversion with Look Ahead.
%
%   `[f,relres,iter,c]=rtisireal(...)` additionally returns an array
%   of residuals `relres`, the number of per-frame iterations done `iter` and the
%   coefficients *c* with the reconstructed phase. The relationship between
%   *f* and *c* is::
%
%     f = idgtreal(c,gd,a,M)
%
%   where *gd* is the canonical dual window obtained by |gabdual|.
%
%   `rtisireal` takes the following addition arguments:
%
%   Algorithm parameters:
%
%     'lookahead',lookahead  Number of lookahead frames. The default value
%                            is `ceil(M/a)-1`.
%
%   !!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!
%   This function requires LTFAT 2.1.2 and above.
%   !!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!
%
%   See also:  idgtreal, dgtreal
%
%   References:
%
%   This function implements algorithm published in:
%
%   Zhu X., Beauregard G. T., Wyse L. L.: Real-Time Signal Estimation From
%   Modified Short-Time Fourier Transform Magnitude Spectra,
%   Audio, Speech, and Language Processing, IEEE Transactions on , vol.15,
%   no.5, pp.1645,1653, July 2007
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
[flags,kv,Ls]=ltfatarghelper({'Ls','maxit'},definput,varargin);

complainif_notposint(kv.maxit,'maxit',mfilename);

% Default number of lookahead frames.
% The kernel size in the horizontal direction is 2*kv.lookahead + 1
if isempty(kv.lookahead)
    kv.lookahead = ceil(M/a)-1;
else
    if kv.lookahead < 0 || kv.lookahead > N - 1
        error('%s: lookahead must be in range [0-%d]',upper(mfilename),N-1);
    end
end

abss = abs(s);
norm_s = norm(abss,'fro');

lookback = ceil(M/a) - 1;
% Analysis window (as array)
gnum = gabwin(g,a,M);
% Synthesis window
gd = gabdual(g,a,M,L);
% ... as array
gdnum = gabwin(gd,a,M);
% Create and overlay several scaled dual windows
wins = repmat(gdnum,1,2*lookback+1);
specg2 = overlayframes(wins,a,M);
% ... and get the asymetric window to be used for the newest lookahead
% frame in 2nd and further iterations
specg2 = (specg2(1:M));
wins(:,1) = 0;
% Do the same thing without contribution of the first window
% This is to be used 
specg1 = overlayframes(wins,a,M);
specg1 = (specg1(1:M));

gnum = fftshift(gnum);
gdnum = fftshift(gdnum);

% Buffer initialization
cframes = zeros(M ,lookback + kv.lookahead+1);
sframes = zeros(M2,kv.lookahead+1);

c = zeros(M2,N);

%recframes = zeros(M,N);

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

% Alternative way of reconstruction from recframes
% 
% f = zeros(L,1);
% idxrange = [0:floor(M/2),-ceil(M/2)+1:-1];
% for n=0:N-1
%     idx = mod(n*a + idxrange,L) + 1;
%     f(idx) = f(idx) + recframes(:,n+1);
% end
% c = dgtreal(f,g,a,M,flags.phase);

f = idgtreal(c,gd,a,M,'timeinv');

if nargout>1
    relres = norm(dgtreal(f,g,a,M,'timeinv')-abss,'fro')/norm_s;
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






