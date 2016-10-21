function [gana,gd] = comp_gsrtisilawins(g,gd,a,M,lookahead)

lookback = ceil(numel(g)/a) -1;

% ... as array
gd = gabwin(gd,a,M);
% Create and overlay several scaled dual windows
wins = repmat(g.*gd*M,1,lookback+1+lookahead);
winsum = comp_overlayframes(wins,a,M);
rellim = 1e-3;
idx = abs(winsum)<rellim;

winsum(idx & winsum > 0) = rellim;
winsum(idx & winsum < 0) = -rellim;
winsum(winsum == 0) = 1;

gana = zeros(numel(g),lookahead+1);
for la = 1:lookahead + 1
    idx = la + lookback - 1;
    gana(:,la) = fftshift(g)./winsum(a*idx+1:a*idx + M);
end

gd = fftshift(gd);
% Reverse order of the windows such that gana(:,1) is analysis window for
% the newest lookahead frame.
gana = gana(:,end:-1:1); 


