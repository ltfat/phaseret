function plotdgtrealphasediff(phase1,phase2,s,thr,a,M,varargin)
%PLOTDGTREALPHASEDIFF Plot DGTREAL phase difference
%   Usage: plotdgtrealphasediff(phase1,phase2,s,thr,a,M)
%
%   Input parameters:
%         phase1  : First phase.  
%         phase2  : Second phase.
%         s       : Magnitude.
%         thr     : Relative threshold
%         a       : Hop size.
%         M       : Number of channels.
%
%   `plotdgtrealphasediff(phase1,phase2,s,thr,a,M)` plots phase difference
%   of *phase1* and *phase2* for lattice given by *a* and *M* taking the 
%   phase wrapping into account. The difference is set to zero for 
%   coefficients with relative magnitude below *thr* i.e. with indices
%   'find(abs(s)<thr*max(abs(s(:))))'. The phase difference is divided by 
%   pi.
%
%   Further, the function accepts any flags and key-value pairs from 
%   |plotdgtreal|. Positional arguments are not recognized (they can 
%   however be passed as key-value arguments).
%

% AUTHOR: Zdenek Prusa

if ~isequal(size(s),size(phase1),size(phase2))
    error('%s: phase1, phase2 and s must have equal dimensions.',...
          upper(mfilename));
end

if ~isscalar(thr)
    error('%s: thr must be a scalar value.',upper(mfilename))
end

d1=abs(modcent(phase1,2*pi)-modcent(phase2,2*pi));
d2=2*pi-d1;
d=min(d1,d2);
d(abs(s)<thr*max(abs(s(:)))) = 0;


d = d/pi;
plotdgtreal(d,a,M,'clim',[0,1],'linabs',varargin{:});



