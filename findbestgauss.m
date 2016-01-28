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

if isempty(atheightrange)
    atheightrange = 0.01:0.001:0.8;
end

w = winwidthatheight(gnum, atheightrange);

L = 10*numel(gnum);

gnum = fir2long(gnum,L);
norms = zeros(size(atheightrange));
tfrs = zeros(size(atheightrange));
for ii=1:numel(atheightrange)
    [gausstmp,tfrs(ii)] = pgauss(L,'inf','width',w(ii),'atheight',atheightrange(ii));
    norms(ii) = norm(gnum-gnum(1)*gausstmp);
end

[~,idx]=min(norms);

atheight = atheightrange(idx);
minorm = norms(idx);
