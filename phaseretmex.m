function phaseretmex(varargin)
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
is_mingw = 0;

if ~exist('comp_rtisilaupdate','file')
  disp(' ');
  disp('--- PHASERET - Phase ReTrieval toolbox. ---');
  disp(' ')
  disp('To start the toolbox, call PHASERETSTART as the first command.');
  disp(' ');
  return;
end;

definput.flags.target = {'auto','lib','mex'};
definput.flags.comptarget = {'release','debug'};
definput.flags.command = {'compile','clean'};
definput.flags.verbosity = {'quiet','verbose'};
flags = ltfatarghelper({},definput,varargin);

do_compilelib = flags.do_lib || flags.do_auto;
do_compilemex = flags.do_mex || flags.do_auto;

if ispc && ~isoctave()
    makecmd = 'mingw32-make';
    is_mingw = 1;
end

try
    if flags.do_compile
        if do_compilelib
            cd([thisdir,filesep,'libltfat']);
            disp('********* Compiling libltfat **********');
            params = ' static NOBLASLAPACK=1';
            if is_mingw, params = [params, ' MINGW=1']; end
            if flags.do_debug, params = [params, ' COMPTARGET=debug']; end
            
            [status,res] = system([makecmd,params]);
            if status ~=0, 
                error(res);
            elseif flags.do_verbose && ~isoctave()
                disp(res);
            end

            cd([thisdir,filesep,'libphaseret']);
            disp('********* Compiling libphaseret **********');
            params = ' static NOBLASLAPACK=1';
            if is_mingw, params = [params, ' MINGW=1']; end
            if flags.do_debug, params = [params, ' COMPTARGET=debug']; end

            [status,res] = system([makecmd,params]);
            if status ~=0
                error(res);
            elseif flags.do_verbose && ~isoctave()
                disp(res);
            end
        end

        if do_compilemex
            cd([thisdir,filesep,'mex']);
            disp('********* Compiling MEX files **********');
            if ~isoctave()
                if ispc
                    params = [' -f Makefile_mingw matlab',...
                    ' MATLABROOT="',matlabroot,'" ARCH=',computer('arch'),...
                    ' EXT=',mexext];
                else
                    params = ' matlab';
                end
            else
                params = ' octave';
            end
            
            if flags.do_debug, params = [params, ' COMPTARGET=debug']; end
            [status,res] = system([makecmd, params]);
            if status ~=0
                error(res);
            elseif flags.do_verbose && ~isoctave()
                disp(res);        
            end
            
        end
    end

    if flags.do_clean
        if do_compilemex
            disp('********* Cleaning MEX files **********');
            cd([thisdir,filesep,'mex']);
            if is_mingw
                [status,res] = system([makecmd,' -f Makefile_mingw clean']);
            else
                [status,res] = system([makecmd,' clean']);
            end
            if status ~=0
                error(res);
            elseif flags.do_verbose && ~isoctave()
                disp(res);
            end
        end

        if do_compilelib
            disp('********* Cleaning libs **********');
            cd([thisdir,filesep,'libphaseret']);
            [status,res] = system([makecmd,' clean']);
            if status ~=0, error(res); end

            cd([thisdir,filesep,'libltfat']);
            [status,res] = system([makecmd,' clean']);
            if status ~=0 
                error(res);
            elseif flags.do_verbose && ~isoctave()
                disp(res);
            end
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

