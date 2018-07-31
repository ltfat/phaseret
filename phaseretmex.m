function phaseretmex(varargin)
%PHASERETMEX  Compiles the MEX interfaces
%   Usage: phaseretmex;
%
%   `phaseretmex` compiles the MEX files. If you have downloaded the binary
%   package, the mex files are already compiled.
%
%   The action of `phaseretmex` is determined by one of the following flags:
%
%     'compile'  Compile stuff. This is the default.
%
%     'clean'    Removes the compiled functions.
%
%   The target to work on is determined by on of the following flags.
%
%   General commands:
%
%     'lib'      Perform action on the LTFAT C library.
%
%     'mex'      Perform action on the mex / oct interfaces.
%
%     'auto'     Choose automatically which targets to work on from the
%                previous ones based on the operation system etc. This is
%                the default.
%
%   Other:
%
%      'verbose' Print action details.
%
%      'debug'   Build a debug version. This will disable compiler
%                optimizations and include debug symbols.

%AUTHOR: Zdenek Prusa

currdir = pwd;
thisdir = fileparts(which(mfilename));
makecmd = 'make';
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

if ispc,   shared_ext = 'dll'; end
if isunix, shared_ext = 'so'; end
if ismac,  shared_ext = 'dylib'; end

try
    if flags.do_compile
        if do_compilelib
            cd([thisdir,filesep,'libltfat']);
            disp('********* Compiling libltfat **********');
            params = ' build/libltfatd.a MODULE=libltfat';
            params = [params, ' NOBLASLAPACK=1 OPTLFLAGS=-DLTFAT_LARGEARRAYS'];
            if is_mingw,
                params = [params, ' MINGW=1 MAKECMDGOALS=static'];
            end
            if flags.do_debug
                params = [params, ' COMPTARGET=debug'];
            end

            [status,res] = system([makecmd,params]);
            if status ~=0,
                error(res);
            elseif flags.do_verbose && ~isoctave()
                disp(res);
            end

            %cd([thisdir,filesep,'libphaseret']);
            disp('********* Compiling libphaseret **********');
            params = ' build/libphaseretd.a MODULE=libphaseret';
            params = [params,' NOBLASLAPACK=1 OPTLFLAGS=-DLTFAT_LARGEARRAYS'];
            if is_mingw
                params = [params, ' MINGW=1 MAKECMDGOALS=static'];
            end
            if flags.do_debug
                params = [params, ' COMPTARGET=debug'];
            end

            [status,res] = system([makecmd,params]);
            resolveres(status,res,flags);
        end

        if do_compilemex
            paramstail = '';
            cd([thisdir,filesep,'mex']);
            disp('********* Compiling MEX files **********');
            if ~isoctave()
                libfftwd = searchforlibfftw(shared_ext);
                paramstail = [' MATLABROOT="',matlabroot,'"',...
                              ' ARCH=',computer('arch'),...
                              ' EXT=',mexext,...
                              ' FFTW_LIBS=',libfftwd];
                
                if ispc
                    params = ' -f Makefile_mingw matlab';
                else
                    [status,res] = system('mex -help');
                    if status ~= 0
                        error('The mex executable in not in the PATH');
                        error('Add %s to your PATH.',[matlabroot,filesep,'bin']);
                    end
                    if ~strcmp(res(1:3),'MEX')
                        [~,res]=system('which mex');
                        error('%s is not a Matlab mex executable!!');
                    end

                    params = ' matlab';
                end
            else
                params = ' octave';
            end

            if flags.do_debug
                params = [params, ' COMPTARGET=debug'];
            end

            if ~isoctave
                matlabversion = version('-release');
                if str2double(matlabversion(1:4)) >= 2018
                    params = [params,' POST2018a=1']; 
                end
            end

            [status,res] = system([makecmd, params, paramstail]);
            resolveres(status,res,flags);
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

            resolveres(status,res,flags);
        end

        if do_compilelib
            disp('********* Cleaning libs **********');
            cd([thisdir,filesep,'libltfat']);
            [status,res] = system([makecmd,' clean']);
            resolveres(status,res,flags);
        end
    end
catch
    % Fallback
    cd(currdir);
    err = lasterror;
    error('Make failed with: \n %s',err.message);
end


cd(currdir);

function libfftw3=searchforlibfftw(sharedExt)
libfftw3base = 'fftw3';
if ispc 
    bp=mfilename('fullpath');
    bp=bp(1:end-length(mfilename));    
    % Search the ltfat/mex lib
    L = dir([bp,'mex',filesep,'*',libfftw3base,'*.',sharedExt]);
    if isempty(L)
        error(['%s: %s could not be found in phaseret/mex subdir.',...
               ' Please download the FFTW dlls and install them.'],...
              upper(mfilename),libfftw3base);
    end
    libfftw3 = L(1).name;
    fprintf('   ...using %s from phaseret/mex.\n',L(1).name);

elseif isunix
     L = dir([matlabroot,filesep,'bin',filesep,computer('arch'),...
              filesep,'*',libfftw3base,'*.',sharedExt,'*']); 

     if isempty(L)
         error('%s: Matlab FFTW libs were not found. Strange.',...
              upper(mfilename));
     end

     libfftw3 = ['"',L(1).folder, filesep, L(1).name,'"'];

     fprintf('   ...using %s from Matlab installation.\n',...
             libfftw3);
end


function isoct=isoctave()
isoct = exist('OCTAVE_VERSION','builtin') ~= 0;

function resolveres(status,res,flags)
if status ~=0
    error(res);
elseif flags.do_verbose && ~isoctave()
    disp(res);
end

