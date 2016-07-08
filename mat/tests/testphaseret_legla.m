function test_failed=test_legla
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
    
    tic;
    chat=legla(s,g,a,M,phaseconv);
    toc;
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('LEGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,relres,iter,f]=legla(s,g,a,M,phaseconv,'modtrunc');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('LEGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,relres,iter,f]=legla(s,g,a,M,phaseconv,'onthefly');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('LEGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,relres,iter,f]=legla(s,g,a,M,phaseconv,'onthefly','modtrunc');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('LEGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,relres,iter,f]=legla(s,g,a,M,phaseconv,'flegla');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FLEGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);   
    
    [chat,relres,iter,f]=legla(s,g,a,M,phaseconv,'flegla','modtrunc');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FLEGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);  
    
    [chat,relres,iter,f]=legla(s,g,a,M,phaseconv,'flegla','onthefly');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FLEGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail); 
    
    [chat,relres,iter,f]=legla(s,g,a,M,phaseconv,'flegla','onthefly','modtrunc');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FLEGLA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail); 
end



function el = getat(collection,id)
if iscell(collection)
    el = collection{id};
else    
    el = collection(id);
end

