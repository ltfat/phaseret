function test_failed=testphaseret_phaseretdemos
%TEST_DEMOS  Test if all the demos runs without errors.

test_failed=0;
matpath = fileparts(which('phaseretstart'));
s=dir([matpath,filesep,'demos',filesep,'demo_*.m']);

for ii=1:numel(s)
     filename = s(ii).name;

     disp(filename);

     % The demo is run in separate function to avoid 
     % variable name clash
     test_failed = test_failed + rundemo(filename(1:end-2));
end


function test_failed = rundemo(demoname)
test_failed = 0;
close all;
try
   eval(demoname);
catch
   test_failed = 1;
   err = lasterror;
   fprintf('DEMO %s FAILED\n with the following reason:',demoname);
   warning(err.message);
   close all;
end
close all;
