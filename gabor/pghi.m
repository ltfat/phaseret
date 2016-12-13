function [c,newphase,usedmask,tgrad,fgrad]=pghi(s,gamma,a,M,varargin)
%PGHI Phase Gradient Heap Integration
%   Usage:  c=pghi(s,gamma,a,M);
%           c=pghi(c,gamma,a,M,mask);
%           c=pghi(c,gamma,a,M,mask,usephase);
%           [c,newphase,usedmask,tgrad,fgrad] = pghi(...);
%
%   Input parameters:
%         s        : Initial coefficients.
%         gamma    : Window "width" factor.
%         a        : Hop factor.
%         M        : Number of channels.
%         mask     : Mask for selecting known phase.
%         usephase : Explicit known phase.
%   Output parameters:
%         c        : Coefficients with the constructed phase.
%         newphase : Just the (unwrapped) phase.
%         usedmask : Mask used for selecting coefficients with the new phase.
%         tgrad    : Relative time phase derivative.
%         fgrad    : Relative frequency phase derivative.
% 
%   `pghi(s,gamma,a,M)` creates complex DGT coefficients from their
%   absolute values *s* using the Phase Gradient Heap Integration algorithm.
%   *s* must have been obtained as::
%
%       c = dgtreal(f,g,a,M,'timeinv');
%       s = abs(c);
%
%   and the algorithm attempts to recover *c*. Parameter *gamma* is window 
%   *g* specific and it can be computed using |pghi_findgamma|.
%
%   `pghi(c,gamma,a,M,mask)` does the same as above except it reuses phase
%   of coefficients *c* for which the corresponding element from *mask*
%   is nonzero.
%
%   `pghi(c,gamma,a,M,mask,usephase)` does the same as above but the known
%   phase *usephase* is passed explicitly.
%
%   This is just a wrapper around |constructphasereal| from LTFAT. Please
%   see its help for more details.
%
%   See also: dgtreal, idgtreal, rtpghi, pghi_findgamma
%
%   References: ltfatnote040 ltfatnote043
%

% AUTHOR: Zdenek Prusa

if ~isscalar(gamma)
    error('%s: gamma must be scalar.',upper(mfilename));
end

N = size(s,2);
L = N*a;
g = {'gauss',gamma/L};

definput.keyvals.tol=[1e-1,1e-10];
definput.keyvals.mask=[];
definput.keyvals.usephase=[];
definput.flags.phase={'timeinv','freqinv'};
[flags,kv]=ltfatarghelper({'mask','usephase'},definput,varargin);

[c,newphase,usedmask,tgrad,fgrad]=...
    constructphasereal(s,g,a,M,'tol',kv.tol,'mask',kv.mask,'usephase',kv.usephase,flags.phase);
