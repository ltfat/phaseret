%DEMO_PGHIOVERLAPS Performance of the PGHI algorithm for several window overlaps 
%
%   The demo performs several reconstructions from the magnitude of the DGT
%   using increasing window overlap: 50%, 75%, 87.5% and 93.75%
%   The number in dB is the spectral convergence of the reconstructed
%   coefficients.
%
%   References: ltfatnote040

% AUTHOR: Zdenek Prusa
%

[f,fs] = cocktailparty;

M = 2048;
g = {'truncgauss',M};
gamma = pghi_findgamma(g); % gamma=lambda*L
stft_magnitude = @(f,a) abs(dgtreal(f,g,a,M,'timeinv'));
inverse_stft = @(c,a) idgtreal(c,{'dual',g},a,M,'timeinv');

% 50% overlap
a = 1024;
s = stft_magnitude(f,a);
c = pghi(s,gamma,a,M);
fhat50p = inverse_stft(c,a);
C50p = magnitudeerrdb(s,stft_magnitude(fhat50p,a));

% 75% overlap
a = 512;
s = stft_magnitude(f,a);
c = pghi(s,gamma,a,M);
fhat75p = inverse_stft(c,a);
C75p = magnitudeerrdb(s,stft_magnitude(fhat75p,a));

% 87.5% overlap
a = 256;
s = stft_magnitude(f,a);
c = pghi(s,gamma,a,M);
fhat875p = inverse_stft(c,a);
C875p = magnitudeerrdb(s,stft_magnitude(fhat875p,a));

% 93.75% overlap
a = 128;
s = stft_magnitude(f,a);
c = pghi(s,gamma,a,M);
fhat9375p = inverse_stft(c,a);
C9375p = magnitudeerrdb(s,stft_magnitude(fhat9375p,a));

disp('To play the original run: soundsc(f,fs);');
fprintf('To play the 1/2 window overlap reconstruction (C=%.2fdB): soundsc(fhat50p,fs);\n',C50p);
fprintf('... 1/4  window overlap (C=%.2fdB): soundsc(fhat75p,fs);\n',C75p);
fprintf('... 1/8  window overlap (C=%.2fdB): soundsc(fhat875p,fs);\n',C875p);
fprintf('... 1/16 window overlap (C=%.2fdB): soundsc(fhat9375p,fs);\n',C9375p);
