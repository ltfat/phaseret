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

if nargin>0
    if strcmpi(target,'clean')
        do_clean = 1;
        do_compilelib = 0;
        do_compilemex = 0;
    end
end

try
    cd([thisdir,filesep,'..']);
    
    if do_compilelib
        disp('********* Compiling library **********');
        [status,res] = system('make lib');
        if status ~=0
            error(res);
        end
    end
    
    if do_compilemex
        disp('********* Compiling MEX files **********');
        if exist('OCTAVE_VERSION','builtin') == 0
            [status,res] = system('make matlab');
        else
            [status,res] = system('make octave');
        end
        if status ~=0
            error(res);
        end
    end
    
    if do_clean
        disp('********* Cleaning MEX files **********');
        [status,res] = system('make clean');
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