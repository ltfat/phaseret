function test_failed=testphaseret_decolbfgs
%TEST_CONSTRUCTPHASE  
%
test_failed = 0;

f = greasy;
a = 128;
M = 1024;
g = {'blackman',M};

for pcId = 1:2
    % Complex case
    phaseconv = getat({'timeinv','freqinv'},pcId);
    
    tra = @(f) dgtreal(f,g,a,M,phaseconv);
    itra = @(c) idgtreal(c,{'dual',g},a,M,phaseconv);
    proj = @(c) tra(itra(c));
    c = tra(f);
    s = abs(c);
    
    
    [chat]=decolbfgs(s,g,a,M,phaseconv,'p',2);
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('DECOLBFGS p=2 %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat]=decolbfgs(s,g,a,M,phaseconv);
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('DECOLBFGS p=2/3 %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
end



function el = getat(collection,id)
if iscell(collection)
    el = collection{id};
else    
    el = collection(id);
end




