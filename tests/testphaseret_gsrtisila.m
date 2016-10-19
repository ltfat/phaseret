function test_failed=testphaseret_gsrtisila
%TEST_CONSTRUCTPHASE
%
test_failed = 0;

f = cocktailparty;
a = 64;
M = 1024;
gl = 512;
h = 0.01;
g = gabwin({'gauss','width',gl,'atheight',h,'inf'},a,M,10*gl);
g = long2fir(g,gl);
gamma = -pi/4*(gl)^2/log(h);

for pcId = 1:2
    % Complex case
    phaseconv = getat({'timeinv','freqinv'},pcId);

    tra = @(f) dgtreal(f,g,a,M,phaseconv);
    itra = @(c) idgtreal(c,{'dual',g},a,M,phaseconv);
    proj = @(c) tra(itra(c));
    c = tra(f);
    s = abs(c);

    [chat]=gsrtisila(s,g,a,M,phaseconv,'maxit',2,'lookahead',3);

    E = magnitudeerrdb(s,proj(chat));
    fail = '';
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('GSRTISILA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);

    [chat]=gsrtisila(s,g,a,M,phaseconv,'unwrap','maxit',2,'lookahead',3);

    E = magnitudeerrdb(s,proj(chat));
    fail = '';
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('GSRTISILA %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);

    [chat]=gsrtisila(s,g,a,M,phaseconv,'spsi','maxit',2,'lookahead',3);

    E = magnitudeerrdb(s,proj(chat));
    fail = '';
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('GSRTISILA lookahead=2 %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);
    
    [chat]=gsrtisila(s,g,a,M,phaseconv,'rtpghi',gamma,'maxit',2,'lookahead',3);

    E = magnitudeerrdb(s,proj(chat));
    fail = '';
    if E>-20
        test_failed = test_failed + 1;
        fail = 'FAILED';
    end

    fprintf('GSRTISILA lookahead=2 %s W=%d E=%.2f %s\n',phaseconv,1,E,fail);

end



function el = getat(collection,id)
if iscell(collection)
    el = collection{id};
else
    el = collection(id);
end




