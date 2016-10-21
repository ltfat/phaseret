function [g,gd,specg1,specg2] = comp_rtisilawins(g,gd,a,M)

% TODO: This only works foe symmetric windows. I was lazy.

lookback = ceil(numel(g)/a) -1;

gd = gabwin(gd,a,M);
% Create and overlay several scaled dual windows
wins = repmat(gd,1,2*lookback+1);
specg2 = comp_overlayframes(wins,a,M);
% ... and get the asymetric window to be used for the newest lookahead
% frame in 2nd and further iterations
specg2 = (specg2(1:M));
wins(:,1) = 0;
% Do the same thing without contribution of the first window
% This is to be used 
specg1 = comp_overlayframes(wins,a,M);
specg1 = (specg1(1:M));

g = fftshift(g);
gd = fftshift(gd);


