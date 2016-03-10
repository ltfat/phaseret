%PHASERET - Phase Retrieval Toolbox
%   Copyright Zdenek Prusa, 2016 
%
%   All algorithms assume the DGT to be computed using the following formula:
%
%   ..              L-1
%   .. c(m+1,n+1) = sum f(l+1)*conj(g(l-a*n+1))*exp(-2*pi*i*m*(l-a*n)/M), 
%   ..              l=0
%
%   .. math:: c\left(m+1,n+1\right)=\sum_{l=0}^{L-1}f(l+1)\overline{g(l-an+1)}e^{-2\pi \mi(l-an+1)m/M}
%
%   where $f$ is a signal of length $L$, $g$ is a window, $M$ is the number of
%   FFT channels, $a$ is time hop factor and $m=0,\ldots,floor(M/2)$ and 
%   $n=0,\ldots,N-1$ and $l-an$ is computed modulo *L*. 
%   Such coefficients can be computed and the signal can be reconstructed 
%   using |dgtreal| and |idgtreal|::
%
%               [c,Ls]=dgtreal(f,g,a,M,'timeinv');
%               [fr]=idgtreal(c,gd,a,M,Ls,'timeinv');
%
%   where $gd$ is canonical dual window obtained by |gabdual|.
%   The phase reconstruction simulation would look like this::
%
%               [c,Ls]=dgtreal(f,g,a,M);
%               [s] = abs(c);
%               [crec] = phaserecalg(...);
%               [fr]=idgtreal(crec,gd,a,M,Ls,'timeinv');
%
%  Implementation of offline algorithms
%    gla     - Griffin-Lim Algorithm
%    legla   - Le Rouxs Griffin-Lim algorithm Modifications
%    pghi        - Phase Gradient Heap Integration
%    decolbfgs   - Decorsiere's phase reconstruction algorithm
%
%  Implementation of online algorithms
%    spsi     - Single Pass Spectrogram Inversion (SPSI)
%    rtisila    - Real-Time Iterative Spectrogram Inversion (RTISI and RTISI-LA)
%    lertisila  - |rtisila| using Le Roux's phase updates
%    rtpghi       - Real-Time Phase Gradient Integration
%
%  Auxilary functions
%    winwidthatheight   - Get window width at specified height
%    findbestgauss      - Find parameters of a Gaussian window closest to a given one
%    findwindowconstant - Find Cg and gamma parameters for |pghi| and |rtpghi|
%    magnitudeerr       - Computes spectral convergence
%    magnitudeerrdb     - |magnitudeerr| in dB
%
%  Plots
%    plotdgtrealphasediff - Plots phase difference
%
%  Demos
%    demo_windowaliasing  - Window aliasing dependence on Gabor sampling lattice
%    demo_vocoderphaseret - Vocoder-based time stretching and pitch shifting
%
%  Block-processing demos
%    demo_blockproc_phaseret    - Comparison of SPSI, RTPGHI and RTISI-LA algorithms, reconstruction only
%    demo_blockproc_phaseretmix - Comparison of SPSI, RTPGHI and RTISI-LA algorithms, delay line mixing
%
%  Note that LTFAT 2.1.2 or newer is required in order to run functions and
%  scripts.





