function [g,gd,specg1,specg2] = comp_rtisilawins(g,gd,a,M)

lookback = ceil(numel(g)/a) -1;

gd = gabwin(gd,a,M);
% Create and overlay several scaled dual windows
wins = repmat(gd,1,2*lookback+1);
specg2 = overlayframes(wins,a,M);
% ... and get the asymetric window to be used for the newest lookahead
% frame in 2nd and further iterations
specg2 = (specg2(1:M));
wins(:,1) = 0;
% Do the same thing without contribution of the first window
% This is to be used 
specg1 = overlayframes(wins,a,M);
specg1 = (specg1(1:M));

g = fftshift(g);
gd = fftshift(gd);


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