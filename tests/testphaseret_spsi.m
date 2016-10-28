function test_failed=testphaseret_spsi
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
    
    [chat]=spsi(s,a,M,phaseconv);
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-10
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('SPSI %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    
    [M2,N] = size(s);
    mask = ones(size(s));
    mask(:,floor(N/2):end) = 0;
    [chat]=spsi(s,a,M,mask,angle(c),phaseconv);
    
    E = magnitudeerrdb(s,proj(chat));
    fail = '';            
    if E>-10
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('SPSI MASK %s W=%d E=%.2f %s\n',phaseconv,1,E,fail); 
end



function el = getat(collection,id)
if iscell(collection)
    el = collection{id};
else    
    el = collection(id);
end

function f = updatef(fpar,f)

