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
%  Implementation of offline algorithms
%    glareal     - Griffin-Lim Algorithm
%    leglareal   - Le Rouxs Griffin-Lim algorithm Modifications
%    pghi        - Phase Gradient Heap Integration
%    decolbfgs   - Decorsiere's phase reconstruction algorithm
%
%  Implementation of online algorithms
%    spsireal     - Single Pass Spectrogram Inversion (SPSI)
%    rtisireal    - Real-Time Iterative Spectrogram Inversion (RTISI and RTISI-LA)
%    lertisireal  - |rtisireal| using Le Roux's phase updates
%    rtpghi       - Real-Time Phase Gradient Integration
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
%  Demos
%    demo_windowaliasing  - Window aliasing dependence on Gabor sampling lattice
%
%  Block-processing demos
%    demo_blockproc_phaseret    - Comparison of SPSI, RTPGHI and RTISI-LA algorithms, reconstruction only
%    demo_blockproc_phaseretmix - Comparison of SPSI, RTPGHI and RTISI-LA algorithms mixing
%
%  Note that LTFAT 2.1.2 or newer is required in order to run functions and
%  scripts.





