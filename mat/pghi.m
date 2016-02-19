function [c,newphase,usedmask,tgrad,fgrad]=pghi(s,g,a,M,varargin)
%PGHI Phase Gradient Heap Integration
%   Usage:  c=pghi(s,g,a,M);
%           c=pghi(s,g,a,M,tol);
%           c=pghi(c,g,a,M,tol,mask);
%           c=pghi(c,g,a,M,tol,mask,usephase);
%           [c,newphase,usedmask,tgrad,fgrad] = pghi(...);
%
%   Input parameters:
%         s        : Initial coefficients.
%         g        : Analysis Gabor window.
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
%   This is just a wrapper around |constructphasereal| from LTFAT. Please
%   see its help for more details. 
%
%   References: ltfatnote040
%

% AUTHOR: Zdenek Prusa

[c,newphase,usedmask,tgrad,fgrad]=constructphasereal(s,g,a,M,varargin{:});