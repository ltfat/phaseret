function chat  = spsi(s,a,M,varargin)
%SPSI Single Pass Spectrogram Inversion (SPSI) for real signals
%   Usage:  c=spsi(s,a,M);
%           c=spsi(s,a,M,mask,phase);   
%
%   Input parameters:
%         s     : $M2 \times N$ array of modulus of Gabor coefficients.
%         a     : Length of time shift.
%         M     : Number of channels.
%   Output parameters:
%         chat   : $M2 \times N$ array of coefficients with reconstructed phase.
%
%   `c=spsi(s,a,M)` returns array of coefficients *c* with reconstructed 
%   phase that the they can be directly used in |idgtreal| such as::
%
%       fhat = idgtreal(...,'timeinv').
%
%   This code was downloaded from
%   http://anclab.org/software/phaserecon/m-files.zip (on 16.9.2015). 
%   Original name of the function was spsi.m
%   It has been modified to work on the Gabor coefficients directly and 
%   +/- pi alternation is no longer necessary.
%
%   References: be15 ltfatnote040
%

%   AUTHORS: Beauregard, G., Harish, M. and Wyse, L.
%   MODIFIED BY: Zdenek Prusa


complainif_notposint(a,'a',mfilename);
complainif_notposint(M,'M',mfilename);

if ~isnumeric(s) || isempty(s)
    error('%s: c must be a numeric array of coefficients.',upper(mfilename));
end

if size(s,3)>1
   error('%s: c cannot be 3dimensional.',upper(mfilename)); 
end

definput.flags.phase={'timeinv','freqinv'};
definput.keyvals.mask = [];
definput.keyvals.phase = [];
[flags,kv,mask,phase]=ltfatarghelper({'mask','phase'},definput,varargin);

M2 = size(s,1); 
M2user = floor(M/2) + 1;

if M2~=M2user
    error('%s: Size of s does not comply with M.',upper(mfilename));
end

if ~isempty(mask)
    if isempty(phase)
        error('%s: mask and phase must be both defined.',upper(mfilename));
    end

    if ~all(cellfun(@(el) isequal(size(el),size(s)),{mask,phase}))
        error('%s: Dimensions of c, mask and phase must be equal.',upper(mfilename));
    end

    % Convert to logical
    % Sanitize mask (anything other than 0 is true)
    mask = cast(mask,'double');
    mask(mask~=0) = 1;
    
    if ~flags.do_timeinv
        % Convert to frequency invariant phase
        phase = angle(phaselockreal(exp(1i*phase),a,M));
    end
end

if isempty(mask)
    chat = comp_spsi(s,a,M);
else
    chat = comp_maskedspsi(s,a,M,mask,phase);
end

if ~flags.do_timeinv
    chat = phaseunlockreal(chat,a,M);
end


