function E = magnitudeerr(target,reconstructed)
%MAGNITUDEERR Computes spectral convergence
%   Usage: E = magnitudeerr(target,reconstructed)
%
%   `magnitudeerr(target,reconstructed)` computes spectral convergence
%   of 2 magnitude spectrograms.
%
%   References: griflim84

E = norm(abs(target)-abs(reconstructed),'fro')/norm(abs(target),'fro');