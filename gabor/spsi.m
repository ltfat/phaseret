function c  = spsi(s,a,M,varargin)
%SPSI Single Pass Spectrogram Inversion (SPSI)
%   Usage:  c=spsi(s,a,M);
%           c=spsi(s,a,M,mask);   
%           c=spsi(s,a,M,mask,usephase);   
%
%   Input parameters:
%         s        : Modulus of Gabor coefficients.
%         a        : Length of time shift.
%         M        : Number of channels.
%         mask     : Mask for selecting known phase.
%         usephase : Explicit known phase.
%   Output parameters:
%         c        : Reconstructed Gabor coefficients.
%
%   `spsi(s,a,M)` creates complex DGT coefficients from their absolute
%   values *s* using the Single Pass Spectrogram Inversion algorithm.
%   *s* must have been obtained as::
%
%       c = dgtreal(f,g,a,M,'timeinv');
%       s = abs(c);
%
%   and the algorithm attempts to recover *c*. Coefficients *c*
%   can be directly used in |idgtreal| such as::
%
%       fhat = idgtreal(...,'timeinv');
%
%   `spsi(c,a,M,mask)` does the same as above except it reuses phase
%   of coefficients *c* for which the corresponding element from *mask*
%   is nonzero.
%
%   `spsi(s,a,M,mask,usephase)` does the same as above but the known
%   phase *usephase* is passed explicitly.
%
%   The original implementation was downloaded from
%   http://anclab.org/software/phaserecon/m-files.zip (on 16.9.2015). 
%   It has been modified to work on the Gabor coefficients directly and,
%   therefore, +/- pi alternation is no longer necessary.
%
%   References: be15 ltfatnote040
%

%   AUTHORS: Beauregard, G., Harish, M. and Wyse, L.
%   MODIFIED BY: Zdenek Prusa


complainif_notposint(a,'a',mfilename);
complainif_notposint(M,'M',mfilename);

if ~isnumeric(s) || isempty(s)
    error('%s: s must be a numeric array of coefficients.',upper(mfilename));
end

if size(s,3)>1
   error('%s: c cannot be 3dimensional.',upper(mfilename)); 
end

definput.flags.phase={'timeinv','freqinv'};
definput.keyvals.mask = [];
definput.keyvals.usephase = [];
[flags,~,mask,usephase]=ltfatarghelper({'mask','usephase'},definput,varargin);

M2 = size(s,1);
M2user = floor(M/2) + 1;
N = size(s,2);

if M2~=M2user
    error('%s: Size of s does not comply with M.',upper(mfilename));
end

if ~isempty(mask)
    if isempty(usephase)
        usephase = angle(s);
    end
    s = abs(s);

    if ~all(cellfun(@(el) isequal(size(el),size(s)),{mask,usephase}))
        error('%s: Dimensions of c, mask and phase must be equal.',upper(mfilename));
    end
    
    % Convert to logical
    % Sanitize mask (anything other than 0 is true)
    mask = cast(mask,'double');
    mask(mask~=0) = 1;
    
    if ~flags.do_timeinv
        % Change to timeinv since spsi expects timeinv
       usephase = angle(phaselockreal(exp(1i*usephase),a,M));                
    end
    
    c = comp_maskedspsi(s,a,M,mask,usephase);
else
    if ~isreal(s) || any(s(:)<0)
        error('%s: s must be real and non-negative when no mask is used.',upper(mfilename));
    end
    
    c = comp_spsi(s,a,M);
end

if ~flags.do_timeinv
    c = phaseunlockreal(c,a,M);
end


