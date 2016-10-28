function [c,newphase,usedmask,tgrad,fgrad]=pghi(s,gamma,a,M,varargin)
%PGHI Phase Gradient Heap Integration
%   Usage:  c=pghi(s,gamma,a,M);
%           c=pghi(s,gamma,a,M,tol);
%           c=pghi(c,gamma,a,M,tol,mask);
%           c=pghi(c,gamma,a,M,tol,mask,usephase);
%           [c,newphase,usedmask,tgrad,fgrad] = pghi(...);
%
%   Input parameters:
%         s        : Initial coefficients.
%         gamma    : Window width factor.
%         a        : Hop factor.
%         M        : Number of channels.
%         tol      : Relative tolerance.
%         mask     : Mask for selecting known phase.
%         usephase : Explicit known phase.
%   Output parameters:
%         c        : Coefficients with the constructed phase.
%         newphase : Just the (unwrapped) phase.
%         usedmask : Mask for selecting coefficients with the new phase.
%         tgrad    : Relative time phase derivative.
%         fgrad    : Relative frequency phase derivative.
% 
%   `pghi(s,gamma,a,M)` creates complex DGT coefficients from their
%   absolute values *s* using the Phase Gradient Heap Integration algorithm.
%   *s* must have been obtained as::
%
%       c = dgtreal(f,g,a,M);
%       s = abs(c);
%
%   and the algorithm attempts to recover *c*. Parameter *gamma* is window 
%   *g* specific and it can be computed using:
%
%   .. gamma = Cg*gl^2
%
%   .. math:: \gamma = C_g \mathit{gl}^2
%
%   where *gl* is the window length and *Cg* is a window specific constant. 
%   Both *Cg* and *gamma* can be obtained by calling |findwindowconstant|.
%
%   This is just a wrapper around |constructphasereal| from LTFAT. Please
%   see its help for more details.
%
%   See also: dgtreal, idgtreal, rtpghi, findwindowconstant
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
[c,newphase,usedmask,tgrad,fgrad]=constructphasereal(s,g,a,M,'timeinv',varargin{:});
