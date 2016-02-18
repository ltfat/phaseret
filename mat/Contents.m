% PHASERET - Phase Retrieval Toolbox
%   Copyright Zdenek Prusa, 2016 
%
%    
%
%   ..              L-1
%   .. c(m+1,n+1) = sum f(l+1)*conj(g(l-a*n+1))*exp(-2*pi*i*m*l/M), 
%   ..              l=0
%
%   .. math:: c\left(m+1,n+1\right)=\sum_{l=0}^{L-1}f(l+1)\overline{g(l-an+1)}e^{-2\pi \mi lm/M}
%            
%  Implementation of algorithms
%    spsireal    - Single Pass Spectrogram Inversion (SPSI)
%    leglareal   - Le Rouxs Griffin-Lim algorithm Modifications
%    rtisireal   - Real-Time Iterative Spectrogram Inversion (RTISI and RTISI-LA)
%    glareal     - Griffin-Lim Algorithm
%
%
%  Auxilary functions
%    winwidthatheight  - Get window width at specified height
%    findbestgauss     - Find parameters of a Gaussian window closest to a given one
%    magnitudeerr      - Computes spectral convergence
%    magnitudeerrdb    - |magnitudeerr| in dB
%
%  Plots
%    plotdgtrealphasediff - Plots phase difference
%
%  Block-processing demos
%    demo_blockproc_phaseret    - Comparison of SPSI, RTPGHI and RTISI-LA algorithms, reconstruction only
%    demo_blockproc_phaseretmix - Comparison of SPSI, RTPGHI and RTISI-LA algorithms mixing
%
%  Note that LTFAT 2.1.2 or newer is required in order to run functions and
%  scripts.





