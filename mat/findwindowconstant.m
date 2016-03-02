function [Cg, gamma] = findwindowconstant( gnum, varargin)
%FINDWINDOWCONSTANT Find window constant for PGHI and RTPGHI
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
%   References: ltfatnote043

% AUTHOR: Zdenek Prusa

definput.keyvals.atheightrange = [];
definput.keyvals.gl = [];
[~,~,gl,atheightrange]=ltfatarghelper({'gl','atheightrange'},definput,varargin);

if ~isvector(gnum) 
    error('%s: Window must be numeric. See FIRWIN and GABWIN.',upper(mfilename))
end

if ~isvector(gnum) || ~isnumeric(gnum)
    gl = numel(gnum);
end

atheight = findbestgauss( gnum, atheightrange);
w = winwidthatheight(gnum, atheight);
Cg = -pi/4*(w/(gl-1))^2/log(atheight);
gamma = Cg*gl^2;

