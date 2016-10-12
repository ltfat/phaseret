function Edb = magnitudeerrdb(target,reconstructed)
%MAGNITUDEERRDB Computes spectral convergence in dB
%   Usage: Edb = magnitudeerr(target,reconstructed)
%
%   `magnitudeerrdb(target,reconstructed)` computes spectral convergence
%   of 2 magnitude spectrograms in dB.
%
%   References: griflim84

Edb = 20*log10(magnitudeerr(target,reconstructed));