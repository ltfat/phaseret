function test_all_phaseret

testprefix = 'testphaseret_';
tests_todo = { 'gla','legla','spsi','rtpghi','pghi','decolbfgs',...
                'rtisila','lertisila','phaseretdemos' };

total_tests_failed=0;
list_of_failed_tests={};

preccell = {'double'};

for precId = 1:length(preccell)
prec = preccell{precId};
    
for ii=1:length(tests_todo)
    testtorun = [testprefix,tests_todo{ii}];
    fprintf('------Running %s------\n',testtorun);
    test_failed=feval([testprefix,tests_todo{ii}]);
    total_tests_failed=total_tests_failed+test_failed;
    if test_failed>0
        list_of_failed_tests{end+1}=[testprefix,tests_todo{ii},' ',prec];
    end;
end;
end

disp(' ');
if total_tests_failed==0
  disp('ALL TESTS PASSED');
else
  s=sprintf('%i TESTS FAILED',total_tests_failed);
  disp(s);
  disp('The following test scripts contained failed tests');
  for ii=1:length(list_of_failed_tests)
    disp(['   ',list_of_failed_tests{ii}]);
  end;
end;
