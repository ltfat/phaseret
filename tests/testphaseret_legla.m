function test_failed=testphaseret_legla
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
    
    [chat,~,relres]=legla(s,g,a,M,phaseconv);
    
    if strcmp(phaseconv,'timeinv')
        % test if timeinv is default
        [chat2,~,relres2,iter2]=legla(s,g,a,M);
        
        E = norm(chat-chat2,'fro');
        
        fail = '';            
        if E>eps
            test_failed = test_failed + 1;
            fail = 'FAILED';
        end
        
        fprintf('GLA timeinv default W=%d %s\n',1,fail);
    end
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20 || abs(E-20*log10(relres(end)))>1
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('LEGLA              %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,~,relres,iter]=legla(s,g,a,M,phaseconv,'modtrunc');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('LEGLA modtrunc     %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,~,relres,iter]=legla(s,g,a,M,phaseconv,'onthefly');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('LEGLA onthefly     %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,~,relres,iter]=legla(s,g,a,M,phaseconv,'onthefly','modtrunc');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('LEGLA onthefly mt  %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat,~,relres,iter]=legla(s,g,a,M,phaseconv,'flegla');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FLEGLA             %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);   
    
    [chat,~,relres,iter]=legla(s,g,a,M,phaseconv,'flegla','modtrunc');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FLEGLA modtrunc    %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);  
    
    [chat,~,relres,iter]=legla(s,g,a,M,phaseconv,'flegla','onthefly');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FLEGLA onthefly    %s W=%d E=%.2f %s\n',phaseconv,1,E,fail); 
    
    [chat,~,relres,iter]=legla(s,g,a,M,phaseconv,'flegla','onthefly','modtrunc');
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('FLEGLA onthefly mt %s W=%d E=%.2f %s\n',phaseconv,1,E,fail); 
end



function el = getat(collection,id)
if iscell(collection)
    el = collection{id};
else    
    el = collection(id);
end

