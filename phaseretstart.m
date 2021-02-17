function phaseretstart(verbose)
%PHASERETSTART Start the Phase Retrival Toolbox
%   Usage: phaseretstart;
%
%   `phaseretstart` starts the PHASERET toolbox. This must be run
%   before using any of the functions in the toolbox.
%
%   The function checks whether the required version of LTFAT is installed
%   and whether MEX files were compiled.
%
%   A typical startup file could look like::
%
%     addpath('/path/to/my/work/ltfat');
%     addpath('/path/to/my/work/phaseret');
%     ltfatstart;
%     phaseretstart;

%AUTHOR: Zdenek Prusa

if nargin < 1
    verbose = 1;
end

requiredLTFAT = '2.1.3';
pgausspath = which('comp_pgauss');
thispath = fileparts(which(mfilename));
basepath = thispath;

phaseretver = deblank(fileread(fullfile(basepath,'phaseret_version')));

% Check if ltfat is started
if isempty(pgausspath)
    error('%s: LTFAT not found. Please install LTFAT %s or newer.',...
    upper(mfilename),requiredLTFAT);
else
    verstr = deblank(fileread([ltfatbasepath,'ltfat_version']));
    if verstr2num(verstr) < verstr2num(requiredLTFAT)
        error('%s: LTFAT %s or newer is required. Detected %s',...
        upper(mfilename),requiredLTFAT,verstr);
    end
end

addpath(thispath);
addpath(fullfile(thispath,'tests'));
addpath(fullfile(thispath,'comp'));
% MEX must be last so that mex files shadow m files
addpath(fullfile(basepath,'mex'));
addpath(fullfile(basepath,'gabor'));
addpath(fullfile(basepath,'demos'));
addpath(fullfile(basepath,'sigproc'));

if verbose
    fprintf('PHASERET version %s. Copyright 2016 Zdenek Prusa.\n',phaseretver);
end

mexcompiled = exist(['phaseret_dummymex.',mexext],'file') == 3;
if ~mexcompiled
     fprintf('MEX files are not compiled. Consider running phaseretmex.\n');
else
    erroccured = 0;
    callerrmsg = '';
    try
        ret = phaseret_dummymex();
    catch
        erroccured = 1;
        err = lasterror;
        callerrmsg = err.message;
    end

    if erroccured || ~strcmp(ret,'mexok')
       error(['MEX files are broken. Consider running ',...
                '`phaseretmex clean` immediatelly.\n%s'],callerrmsg); 

    end
end


function num = verstr2num(verstr)
vercell = textscan(verstr,'%s','Delimiter','.');
vercell = vercell{1};
vercell = vercell(~cellfun(@isempty,vercell));
vercell = vercell(end:-1:1);
num = 0;
for ii=1:numel(vercell)
     num = num + 1000^ii*str2double(vercell{ii});
end

