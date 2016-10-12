function phaseretmex(target)
%PHASERETMEX  Compiles the MEX interfaces
%   Usage: phaseretmex;
%
%   `phaseretmex` compiles the MEX files. If you have downloaded the binary
%   package, the mex files are already compiled.

currdir = pwd;
thisdir = fileparts(which(mfilename));
do_clean = 0;
do_compilelib = 1;
do_compilemex = 1;
makecmd = 'make';
fftwlibs = '';

if ispc && ~isoctave()
    makecmd = 'mingw32-make';
end

if nargin>0
    if strcmpi(target,'clean')
        do_clean = 1;
        do_compilelib = 0;
        do_compilemex = 0;
    end
end

try

    if do_compilelib
        cd([thisdir,filesep,'libltfat']);
        disp('********* Compiling libltfat **********');
        [status,res] = system([makecmd,' static NOBLASLAPACK=1']);
        if status ~=0
            error(res);
        end
        
        cd([thisdir,filesep,'libphaseret']);
        disp('********* Compiling libphaseret **********');
        [status,res] = system([makecmd,' static']);
        if status ~=0
            error(res);
        end
    end

    if do_compilemex
        cd([thisdir,filesep,'mex']);

        disp('********* Compiling MEX files **********');
        if ~isoctave()
            if ispc
                [status,res] = system([makecmd,' -f Makefile_mingw matlab',...
                ' MATLABROOT="',matlabroot,'" ARCH=',computer('arch'),...
                ' EXT=',mexext]);
            else
                [status,res] = system([makecmd,' matlab']);
            end
        else
            [status,res] = system([makecmd, ' octave']);
        end
        if status ~=0
            error(res);
        end
    end

    if do_clean
        cd([thisdir,filesep,'mex']);
        disp('********* Cleaning MEX files **********');
        if ispc && ~isoctave()
            [status,res] = system([makecmd,' -f Makefile_mingw clean']);
        else
            [status,res] = system([makecmd,' clean']);
        end
        if status ~=0
            error(res);
        end

        cd([thisdir,filesep,'libphaseret']);
        disp('********* Cleaning lib **********');
        [status,res] = system([makecmd,' cleanlib']);
        if status ~=0
            error(res);
        end
    end
catch
    cd(currdir);
    err = lasterror;
    error('Make failed with: \n %s',err.message);
end


cd(currdir);


function isoct=isoctave()
isoct = exist('OCTAVE_VERSION','builtin') ~= 0;
