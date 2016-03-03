function test_failed=test_gla
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
    
    [chat,relres,iter,f]=gla(s,g,a,M,phaseconv);
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('GLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,relres,iter,f]=gla(s,g,a,M,phaseconv,'fgla');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);   
end



function el = getat(collection,id)
if iscell(collection)
    el = collection{id};
else    
    el = collection(id);
end

function f = updatef(fpar,f)

