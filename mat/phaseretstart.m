function phaseretstart(verbose)
%PHASERETSTART Start the Phase Retrival Toolbox
%   Usage: phaseretstart;
%
%   `phaseretstart` starts the PHASERET toolbox. This must be run
%   before using any of the functions in the toolbox.
%
%   The function adds [basepath] and [basepath]/../mex directories to path,
%   where basepath is the directory where phaseretstart resides in.
%   The function also checks whether LTFAT is installed and whether MEX 
%   files were compiled.
%

if nargin < 1
    verbose = 1;
end

requiredLTFAT = '2.1.2';
pgausspath = which('comp_pgauss');
thispath = fileparts(which(mfilename));
basepath = thispath(1:end-4);

phaseretver = deblank(fileread(fullfile(basepath,'phaseret_version')));

% Check if ltfat is started
if isempty(pgausspath)
    error('%s: LTFAT not found. Please install LTFAT 2.1.2 or newer.',...
    upper(mfilename))
else
    verstr = deblank(fileread([ltfatbasepath,'ltfat_version']));
    if verstr2num(verstr) < verstr2num(requiredLTFAT)
        error('%s: LTFAT 2.1.2 or newer is required. Detected %s',...
        upper(mfilename),verstr);
    end
end

addpath(thispath);
addpath(fullfile(basepath,'mex'));
addpath(fullfile(thispath,'tests'));

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

