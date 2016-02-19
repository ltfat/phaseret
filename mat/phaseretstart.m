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

requiredLTFAT = '2.1.1';
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

if verbose
    fprintf('PHASERET version %s. Copyright 2016 Zdenek Prusa.\n',phaseretver);

    mexcompiled = exist(['comp_rtisilaupdate.',mexext],'file') == 3;
    if ~mexcompiled
        fprintf('MEX files are not compiled. Consider running phaseretmex.\n');
    end
end

function num = verstr2num(verstr)
vercell = strsplit(verstr,'.');
vercell = vercell(~cellfun(@isempty,vercell));
vercell = vercell(end:-1:1);
num = 0;
for ii=1:numel(vercell)
     num = num + 100^ii*str2double(vercell{ii});
end

