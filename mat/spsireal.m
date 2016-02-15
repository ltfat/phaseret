function chat  = spsireal(c,a,M,varargin)
%SPSIREAL Single Pass Spectrogram Inversion (SPSI) for real signals
%   Usage:  c=spsireal(s,a,M);
%           c=spsireal(s,a,M,mask);   
%
%   Input parameters:
%         s     : $M2 \times N$ array of modulus of Gabor coefficients.
%         a     : Length of time shift.
%         M     : Number of channels.
%   Output parameters:
%         chat   : $M2 \times N$ array of coefficients with reconstructed phase.
%
%   `c=spsireal(s,a,M)` returns array of coefficients *c* with reconstructed 
%   frequency invariant phase (the default convention in |dgtreal|) such 
%   that the they can be directly used in |idgtreal|.
%
%   `c=spsireal(s,a,M,'timeinv')` returns coefficients with time invariant
%   phase. To reconstruct the signal function |idgtreal| must be called
%   with 'timeinv' i.e::
%   
%       fhat = idgtreal(...,'timeinv').
%
%   This code was downloaded from
%   http://anclab.org/software/phaserecon/m-files.zip (on 16.9.2015). 
%   Original name of the function was spsi.m
%   It has been modified to work on the Gabor coefficients directly and 
%   +/- pi alternation is no longer necessary.
%
%   Examples:
%   ---------
%
%   Reconstructing phase for greasy
%   :::
%   
%   
%   See also:  idgtreal, dgtreal
%
%   References: 
%
%   Beauregard, G., Harish, M. and Wyse, L. (2015), Single Pass Spectrogram
%   Inversion, in Proceedings of the IEEE International Conference on 
%   Digital Signal Processing. Singapore, 2015. 
%
%   Prusa, Z., Balazs, P. and Soendergaard, P. L. (2016), A Non-iterative 
%   Method for STFT Phase (Re)Construction, ...

%   AUTHORS: Beauregard, G., Harish, M. and Wyse, L.
%   MODIFIED BY: Zdenek Prusa


complainif_notposint(a,'a',mfilename);
complainif_notposint(M,'M',mfilename);

if ~isnumeric(c) || isempty(c)
    error('%s: c must be a numeric array of coefficients.',upper(mfilename));
end

if size(c,3)>1
   error('%s: c cannot be 3dimensional.',upper(mfilename)); 
end

definput.flags.phase={'freqinv','timeinv'};
definput.keyvals.mask = [];
[flags,kv]=ltfatarghelper({},definput,varargin);
mask = kv.mask;

if ~isempty(mask)
    % Convert to logical
    mask = kv.mask ~= 0;
end

[M2,N,W] = size(c); 
M2user = floor(M/2) + 1;

if M2~=M2user
    error('%s: Size of s does not comply with M.',upper(mfilename));
end

chat = comp_spsireal(c,a,M,mask);


if ~flags.do_timeinv
    % Convert to frequency invariant phase
    N = size(chat,2);
    L = N*a;
    b = L/M;
    M2=floor(M/2)+1;
    N=size(chat,2);

    TimeInd = (0:(N-1))/N;
    FreqInd = (0:(M2-1))*b;

    phase = FreqInd'*TimeInd;
    phase = exp(-2*1i*pi*phase);

    % Handle multisignals
    chat = bsxfun(@times,chat,phase);
end


