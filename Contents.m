%
%   Note that LTFAT 2.1.2 or newer is required in order to run functions and 
%   scripts.
% 
%   Copyright Zdenek Prusa, 2015
% 
%   Implementation of algorithms (not included in LTFAT)
% 
%      spsireal   - Single Pass Spectrogram Inversion (SPSI)
%      leglareal  - Le Rouxs Griffin-Lim algorithm Modifications
%      rtisireal  - Real-Time Iterative Spectrogram Inversion (RTISI and RTISI-LA)
% 
%   Scripts reproducing tables
% 
%      repr_table1  - Comparison with SPSI, MOCHA-TIMIT database
%      repr_table2  - Comparison with SPSI, EBU SQAM database
% 
%   Scripts reproducing figures
% 
%      repr_fig1    - Phase difference plots for increasing hop factor
%      repr_fig2    - Phase difference plots for 3 different algorithms
%      repr_fig3    - Comparison with iterative algorithms, MOCHA-TIMIT database
%      repr_fig4    - Comparison with iterative algorithms, EBU SQAM database
% 
%   Auxilary functions
% 
%      winwidthatheight  - Get window width at specified height
%      findbestgauss     - Find parameters of a Gaussian window closest to a given one
%      magnitudeerr      - Computes spectral convergence
%      magnitudeerrdb    - |magnitudeerr| in dB
% 
%   Plots
%      plotdgtrealphasediff - Plots phase difference 
% 
% 
%   Demos
%      demo_equalerror   - 
%      demo_modified     -
%      demo_pitchshift   -
%      
    


    
    
