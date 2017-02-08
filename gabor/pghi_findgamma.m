function [gamma, Cg, gl] = pghi_findgamma( g, varargin)
%PGHI_FINDGAMMA Find window constant for PGHI and RTPGHI
%   Usage: gamma = pghi_findgamma({firwinname,gl})
%          gamma = pghi_findgamma(g,a,M)
%          gamma = pghi_findgamma(g,a,M,L)
%          [gamma, Cg] = pghi_findgamma(...)
%
%   Input parameters:
%         gnum     : Window.
%         gl       : Length of the support of the window.
%   Output parameters:
%         gamma    : Parameter for PGHI and RTPGHI
%         Cg       : Window constant
%
%   `pghi_findgamma({firwinname,gl})` returns parameter *gamma*, for which the
%   Gaussian window given as:
%
%   .. g = exp(-pi*l^2/gamma)
%
%   .. math:: g=exp(-\pi\frac{l^2}{\gamma})
%
%   is closest to peak-normalized window *firwinname* from |firwin|.
%   The parameter is precomputed so the search will not be done.
%
%   `[gamma,Cg] = pghi_findgamma(...)` additionaly returns parameter
%   *Cg* which is window constant and is used to compute gamma such as:
%
%   .. gamma = Cg*gl^2
%
%   .. math:: \gamma=C_g \mathit{gl}^2}
%
%   where *gl* is the length of the window support.
%
%   Note that the relationship between *gamma* and *tfr* from |pgauss| is:
%
%   .. tfr = gamma/L
%
%   .. math:: tfr=gamma/L
%
%   where *L* is the DGT length.
%
%   Additional parameters:
%   ----------------------
%
%   'search'             Do the search even for precomputed windows.
%
%   References: ltfatnote043

% AUTHOR: Zdenek Prusa

definput.keyvals.atheightrange = [];
definput.keyvals.a = [];
definput.keyvals.M = [];
definput.keyvals.L = [];
definput.flags.method = {'precomputed','search'};
[flags,kv]=ltfatarghelper({'a','M','L'},definput,varargin);

wins = getfield(arg_firwin,'flags','wintype');

if ischar(g) || (iscell(g) && ischar(g{1})) && ~flags.do_search
    winname = g;
    if iscell(g), winname = g{1}; end;
        switch winname
        case wins
            if iscell(g) && numel(g)>1 && isnumeric(g{2})
                gl = g{2}; 
            elseif ~isempty(kv.M)
                gl = kv.M;
            else
                error('%s: Window length is unspecified.',upper(mfilename));
            end

            try
                [gamma,Cg] = precomputed_gamma(winname,gl);
                return; % Return immediatelly
            catch
            end
        case 'gauss'
            Cg = nan;
            if iscell(g) && numel(g)>1 && isnumeric(g{2})
                if isempty(kv.L)
                    error('%s: Gaussian window requires L to be specified',...
                    upper(mfilename));
                end
                gamma = g{2}*kv.L;
            else
                if any(cellfun(@isempty,{kv.a,kv.M}))
                    error('%s: a and M must be defined for this window.',upper(mfilename));
                end
                gamma = kv.a*kv.M;
            end
            return;
    end
end

if ~isnumeric(g)
    if any(cellfun(@isempty,{kv.a,kv.M}))
        error('%s: a and M must be defined for this window.',upper(mfilename));
    end
    g = gabwin(g,kv.a,kv.M,kv.L);
end

gl = round(winwidthatheight(g, 1e-10));
g = long2fir(g,gl);

atheight = findbestgauss( g, kv.atheightrange);
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
        if isempty(ind2)
            width(ii) = gl;
        else
            rest = 1-(fracofmax-gnum(ind2))/(gnum(ind1)-gnum(ind2));
            width(ii) = 2*(ind1+rest-1);
        end
    else
        width(ii) = 2*(ind-1);
    end
end


function [gamma,Cg] = precomputed_gamma(g,gl)

switch g
case {'hann','hanning','nuttall10'}
        Cg = 0.25645;
    case {'sqrthann','cosine','sine'}
        Cg = 0.41532;
    case {'hamming'}
        Cg = 0.29794;
    case {'nuttall01'}
        Cg = 0.29610;
    case {'tria','triangular','bartlett'}
        Cg = 0.27561;
    case {'sqrttria'}
        Cg = 0.48068;
    case {'blackman'}
        Cg = 0.17954;
    case {'blackman2'}
        Cg = 0.18465;
    case {'nuttall','nuttall12'}
        Cg = 0.12807;
    case {'ogg','itersine'}
        Cg = 0.35744;
    case {'nuttall20'}
        Cg = 0.14315;
    case {'nuttall11'}
        Cg = 0.17001;
    case {'nuttall02'}
        Cg = 0.18284;
    case {'nuttall30'}
        Cg = 0.09895;
    case {'nuttall21'}
        Cg = 0.11636;
    case {'nuttall03'}
        Cg = 0.13369;
    case {'truncgauss'}
        Cg = 0.17054704423023;
    otherwise
        error('Unsupported FIR window type');
end

gamma = Cg*gl^2;



