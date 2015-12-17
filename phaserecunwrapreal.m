function [chat,tonset] = phaserecunwrapreal(c,g,a,M,varargin)
%   Usage:  c=spsireal(s,a,M);
%
%   Input parameters:
%         s     : $M2 \times N$ array of modulus of Gabor coefficients.
%         g     : Window function (in LTFAT format)
%         a     : Length of time shift.
%         M     : Number of channels.
%   Output parameters:
%         chat   : $M2 \times N$ array of coefficients with reconstructed phase.
%         
%         f_centr: Length $N$ cell-array, each element is a vector with 
%                  peak frequencies fom interval [0,M2-1].
%         f_inf  : $M2 \times N$ array of indices of 
%
%   Caution! This code depends on the Signal Processing Toolbox
%
%   This code was downloaded from
%   http://www.perso.telecom-paristech.fr/~magron/ressources/ph_unwrapping_matlab.7z
%   (on 12.11.2015). 
%   Original name of the function was phase_unwrap.m
%   
%   It has been heavily modified
%
%   See also:  idgtreal, dgtreal
%
%   References:
%
%   P. Magron, R. Badeau, B. David. Phase reconstruction of spectrograms 
%   with linear unwrapping: application to audio signal restoration, in 
%   Proc. of the European Signal Processing Conference (EUSIPCO), Nice, 
%   France, 2015.

%   AUTHORS: Paul Magron, Feb 2015
%   MODIFIED BY: Zdenek Prusa


definput.flags.phase={'freqinv','timeinv'};
[flags,kv]=ltfatarghelper({},definput,varargin);

s = abs(c);
[M2,N] = size(s);
L = N*a;
M2user = floor(M/2) + 1;
if M2~=M2user
   error('M does not comply with the size of c.'); 
end

g = gabwin(g,a,M,L);

% Onset detection
tonset = detect_onset_frames(s);
tonset = tonset(:)';

mask = zeros(M2,N); mask = mask~=0;
newphase = zeros(M2,N);

% Reconstruction of vertical frames

% Here we construct frequency invariant phase since it only depends
% on the absolute position of the impuse. 
% The time invariant phase depends on a distance of the window
% from the absolute position of the impulse.

% The impulse position is frequency dependent and it can vary
% slightly.

% Window effective support i.e. width at the 0.01 height
wef = ceil(winwidthatheight(g,0.2));

r = ceil(wef/(2*a));
range = -r:r;

for oId = 1:numel(tonset)

    onsetRange = tonset(oId) + range;
    onsetRange(onsetRange<2) = []; 
    onsetRange(onsetRange>N-1) = []; 

    [~,maxId] = max(s(:,onsetRange),[],2);

    maxId = onsetRange(maxId);
    for m=2:numel(maxId)
        n0frac = qint(20*log10(1+s(m,maxId(m)-1)),...
                      20*log10(1+s(m,maxId(m))),...
                      20*log10(1+s(m,maxId(m)+1)))+maxId(m);
        n0frac = min([n0frac,N]); n0frac = max([n0frac,1]);
        n0range = maxId(m) + range;
        %n0range = onsetRange;
        n0range(n0range<1) = []; n0range(n0range>N) = [];               

        phase_tmp = -2*pi*(n0frac-1)*a/M;

        newphase(m,n0range) = newphase(m-1,round(n0frac)) + phase_tmp;
        mask(m,n0range) = 1;
    end
end

% Freqinv to timeinv
TimeInd = (0:(N-1))*a;
FreqInd = (0:(M2-1))/M;
phase = FreqInd'*TimeInd;
phase = 2*pi*phase;
newphase(mask) = newphase(mask) + phase(mask);


chat = abs(s);
chat(mask) = chat(mask).*exp(1i*newphase(mask));
chat = spsireal(chat, a, M , 'timeinv','mask',mask);

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



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  SUPPORT FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,b,a] = qint(ym1,y0,yp1)
den = 2*(2*y0 - yp1 - ym1);
if den == 0
    p = 0;
    return;
end

p = (yp1 - ym1)./den;
b = y0 - 0.25.*(ym1-yp1).*p;
a = 0.5*(ym1 - 2*y0 + yp1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ONSET DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Onset frames detection from a spectrogram
%
% Inputs :
%     V : spectrogram (F*T matrix)
%     Fs : sample rate
%     w : analysis window
%     hop : hop size
%
% Output :
%     tonset : vector of onset frames

function [tonset] = detect_onset_frames(s)


sdiff = pderiv(s,2);
v0 = sum((sdiff + abs(sdiff))/2);
v0 = v0 - mean(v0);
v0 = v0/std(v0);


[~,ind] = findpeaks(v0,'MINPEAKHEIGHT',max(v0)*0.1);
t0 = ind'; t0 = t0(1:end);
tonset = t0;

end





