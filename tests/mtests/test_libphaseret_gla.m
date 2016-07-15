clear all;
f = greasy;
a = 64;
M = 512;
M2 = floor(M/2) + 1; 
gl = 512;
L = dgtlength(numel(f),a,M);
g = firwin('hann',gl);
gd = long2fir(gabdual(g,a,M),gl);
N = L/a;
maxit = 100;

corig = dgtreal(f,{'hann',gl},a,M,'timeinv');
s = abs(corig);

cout = zeros(2*M2,N);
coutPtr = libpointer('doublePtr',cout);

calllib('libphaseret','gla',s,g,gl,L,1,a,M,maxit,coutPtr);

cout2 = interleaved2complex(coutPtr.Value);

frec = idgtreal(cout2,{'dual',{'hann',gl}},a,M,'timeinv');

s2 = dgtreal(frec,{'hann',gl},a,M,'timeinv');
magnitudeerrdb(s,s2)




