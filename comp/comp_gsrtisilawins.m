function [gana,gd] = comp_gsrtisilawins(g,gd,a,M,lookahead)

lookback = ceil(numel(g)/a) -1;

% ... as array
gd = gabwin(gd,a,M);
% Create and overlay several scaled dual windows
wins = repmat(g.*gd*M,1,lookback+1+lookahead);
winsum = overlayframes(wins,a,M);
rellim = 1e-3;
idx = abs(winsum)<rellim;
winsum(idx & winsum > 0) = rellim;
winsum(idx & winsum == 0) = 1;
winsum(idx & winsum < 0) = -rellim;

gana = zeros(numel(g),lookahead+1);
for la = 1:lookahead + 1
    idx = la + lookback - 1;
    gana(:,la) = fftshift(g)./winsum(a*idx+1:a*idx + M);
end

gd = fftshift(gd);
% Reverse order of the windows such that gana(:,1) is analysis window for
% the newest lookahead frame.
gana = gana(:,end:-1:1); 

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

