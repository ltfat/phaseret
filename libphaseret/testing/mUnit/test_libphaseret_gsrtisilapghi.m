
f = gspi;
a = 128;
M = 1024;
M2 = floor(M/2) + 1; 
gl = 1024;
L = dgtlength(numel(f),a,M);
g = firwin('hann',gl);
gamma = 0.25645*gl^2;
gd = long2fir(gabdual(g,a,M),gl);
N = L/a;
lookahead = 10;
maxit = 10;

corig = dgtreal(f,{'hann',gl},a,M,'timeinv');
s = abs(corig);

cout = zeros(2*M2,N);
coutPtr = libpointer('doublePtr',cout);

calllib('libphaseret','phaseret_gsrtisilapghioffline_d',s,g,L,gl,1,a,M,lookahead,maxit,...
                      gamma, 1e-6, 0, coutPtr);

cout2 = interleaved2complex(coutPtr.Value);

frec = idgtreal(cout2,{'dual',{'hann',gl}},a,M,'timeinv');

s2 = abs(dgtreal(frec,{'hann',gl},a,M,'timeinv'));
magnitudeerrdb(s,s2)



c=gsrtisila(s,g,a,M,'lookahead',lookahead,'maxit',maxit,'timeinv','rtpghi',{gamma});

frec = idgtreal(c,{'dual',{'hann',gl}},a,M,'timeinv');
s3 = abs(dgtreal(frec,{'hann',gl},a,M,'timeinv'));
magnitudeerrdb(s,s3)


