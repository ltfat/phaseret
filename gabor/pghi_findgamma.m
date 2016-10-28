function [gamma, Cg] = pghi_findgamma( g, varargin)
%PGHI_FINDWINDOWCONSTANT Find window constant for PGHI and RTPGHI
%   Usage: Cg = findwindowconstant(gnum)
%          Cg = findwindowconstant(gnum,gl)
%          [Cg, gamma] = findwindowconstant(...)
%
%   Input parameters:
%         gnum     : Window.
%         gl       : Length of the support of the window.
%   Output parameters:
%         Cg       : Window constant
%         gama     : Parameter for PGHI and RTPGHI
%
%   `Cg = findwindowconstant(gnum)` does a heuristic search for the 
%   parameter *Cg*, for which the Gaussian window given as:
%   
%   .. g = exp(-pi*l^2/(Cg*gl^2))
%   
%   .. math:: g=exp(-\pi\frac{l^2}{C_g \mathit{gl}^2})
%
%   is closest to peak-normalized window *gnum*, where *gl* is length 
%   of its support.
%
%   `Cg = findwindowconstant(gnum,gl)` works as before but uses explicitly
%   given *gl*. This is usefull when e.g. *gnum* is zero padded.
%
%   `[Cg,gamma] = findwindowconstant(...)` additionaly returns parameter
%   *gamma* which is equal to:
%   
%   .. gamma = Cg*gl^2
%   
%   .. math:: g=C_g \mathit{gl}^2}
%
%   Note that the relationship between *gamma* and *tfr* from |pgauss| is:
%
%   .. tfr = gamma/L
%   
%   .. math:: tfr=gamma/L
%
%   where *L* is the DGT length.
%
%   References: ltfatnote043

% AUTHOR: Zdenek Prusa

definput.keyvals.atheightrange = [];
definput.keyvals.gl = [];
[~,~,gl,atheightrange]=ltfatarghelper({'gl','atheightrange'},definput,varargin);

if ~isvector(g) 
    error('%s: Window must be numeric. See FIRWIN and GABWIN.',upper(mfilename))
end

if isempty(gl) || ~isvector(g) || ~isnumeric(g)
    gl = numel(g);
end

atheight = findbestgauss( g, atheightrange);
w = winwidthatheight(g, atheight);

Cg = -pi/4*(w/(gl-1))^2/log(atheight);
gamma = Cg*gl^2;


function [atheight,minorm] = findbestgauss( gnum , varargin)
%FINDBESTGAUSS Find Gaussian window closest to the given window
%   Usage: atheight = winwidthatheight(gnum)
%
%   Input parameters:
%         gnum     : Window.
%   Output parameters:
%         atheight : Relative height where both windows have the same width.
%         minorm   : Error of the windows
%
%   `atheight = findbestgauss(gnum)` searches for a Gaussian window
%   which is closest to window *gnum* and returns a relative height *atheight*
%   at which both windows have the same width. *gnum* must be a numeric
%   vector returned from from |gabwin| or |firwin|. The function does a
%   simple heuristic search. Nothing fancy.
%   
%   Examples:
%   ---------
%
%   The following example shows how to use this function to determine
%   parameters of the Gaussian window closest to the Hann window:::
%
%       % Create a probe
%       gnum = firwin('hann',1024);
%       atheight = findbestgauss(gnum)
%       % ...
%       % Elsewhere we want parameters of a Gaussian window for different
%       % Hann window
%       ghann = firwin('hann',2048);
%       width = winwidthatheight(ghann,atheight);
%       L = 4*2048;
%       ggauss = pgauss(L,'width',width,'atheight',atheight);
%       plot(0:L-1,normalize([fir2long(ghann,L),ggauss],'inf'));
%       hold on;
%       lhandle = line([1,L],[atheight,atheight]);
%       set(lhandle,'LineStyle','--'); set(lhandle,'Color','k');
%       hold off;
%
%       % The following can be directly used in |dgtreal| and |constructphasereal|
%       g = {'gauss','width',width,'atheight',atheight};
%             

% AUTHOR: Zdenek Prusa

definput.keyvals.atheightrange = [];
[~,~,atheightrange]=ltfatarghelper({'atheightrange'},definput,varargin);

if ~isvector(gnum) || ~isnumeric(gnum)
    error('%s: Window must be numeric. See FIRWIN and GABWIN.',upper(mfilename))
end

if isempty(atheightrange)
    atheightrange = 0.01:0.001:0.8;
end

w = winwidthatheight(gnum, atheightrange);

L = 10*numel(gnum);

gnum = fir2long(normalize(gnum,'inf'),L);
norms = zeros(size(atheightrange));
tfrs = zeros(size(atheightrange));
for ii=1:numel(atheightrange)
    [gausstmp,tfrs(ii)] = pgauss(L,'inf','width',w(ii),'atheight',atheightrange(ii));
    norms(ii) = norm(gnum-gausstmp);
end

[~,idx]=min(norms);

atheight = atheightrange(idx);
minorm = norms(idx);


function width = winwidthatheight(gnum,atheight)
%WINWIDTHATHEIGHT Window width at height
%   Usage: width = winwidthatheight(gnum, height)
%
%   Input parameters:
%         gnum      : Window.
%         atheight  : Relative height.
%   Output parameters:
%         width   : Window width in samples.
%
%   `winwidthatheight(gnum,atheight)` computes width of a window *gnum* at
%   the relative height *atheight*. *gnum* must be a numeric vector as
%   returned from |gabwin|. If *atheight* is an array, width will have the
%   same shape with correcpondng values.
%

% AUTHOR: Zdenek Prusa

if ~isnumeric(gnum) || isempty(gnum) || ~isvector(gnum)
    error('%s: gnum must be a numeric vector.', upper(mfilename));
end

if isempty(atheight) || any(atheight) > 1 || any(atheight) < 0
    error('%s: h must be in the interval [0-1].', upper(mfilename));
end

width = zeros(size(atheight));
for ii=1:numel(atheight)
    gl = numel(gnum);
    gmax = max(gnum);
    frac=  1/atheight(ii);
    fracofmax = gmax/frac;
    
    
    ind =find(gnum(1:floor(gl/2)+1)==fracofmax,1,'first');
    if isempty(ind)
        %There is no sample exactly half of the height
        ind1 = find(gnum(1:floor(gl/2)+1)>fracofmax,1,'last');
        ind2 = find(gnum(1:floor(gl/2)+1)<fracofmax,1,'first');
        rest = 1-(fracofmax-gnum(ind2))/(gnum(ind1)-gnum(ind2));
        width(ii) = 2*(ind1+rest-1);
    else
        width(ii) = 2*(ind-1);
    end
end


