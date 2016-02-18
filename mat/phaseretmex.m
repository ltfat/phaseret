function phaseretmex(target)
%PHASERETMEX  Compiles the MEX interfaces
%   Usage: phaseretmex;
%
%   `phaseretmex` compiles the MEX files. If you have downloaded the binary
%   package, the mex files are already compiled.

currdir = pwd;
thisdir = fileparts(which(mfilename));

if nargin>0
    if strcmpi(target,'clean')
        
    end
else
    try
        cd([thisdir,filesep,'..']);
        system('make lib');
        if exist('OCTAVE_VERSION','builtin') == 0
            system('make matlab');
        else
            system('make octave');
        end
    catch
        cd(currdir);
        err = lasterror;
        error('Compilation failed with: \n %s',err.message);
    end
end

 cd(currdir);