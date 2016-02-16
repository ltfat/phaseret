%DEMO_EQUALERROR 
%
%   This demo allows comparison of phase reconstruction from 
%   |constructphasereal| and from iterative algorithms. The algorithms
%   are stopped such that the error is equal to the error achieved by
%   |constructphasereal|.
%

[f,fs] = wavload(['Databases',filesep,'SQAM',filesep,'08.wav']);
f = f(1:10*fs,1);
Ls = numel(f);

a = 256;
M = 2048;
gl = M;
thr = [1e-1,1e-10];
thrplot = 1e-4;

L=dgtlength(Ls,a,M);
b = L/M;
N = L/a;

g = 'gauss';
ggauss = g;
[~,info]=gabwin(g,a,M,L);

gamma = 0.9;
graycm = flipud(gray.^gamma);

c = dgtreal(f,g,a,M,'timeinv');
s = abs(c);

chat = constructphasereal(abs(c),ggauss,a,M,thr,'timeinv');

fhatHeapint = idgtreal(chat,{'dual',g},a,M,numel(f),'timeinv');
nextprojc = dgtreal(fhatHeapint,g,a,M,'timeinv');
EheapintDB = magnitudeerrdb(s,nextprojc)
Eheapint = magnitudeerr(s,nextprojc);

figure(1);clf;
plotdgtrealphasediff(angle(chat),angle(c),s,thrplot,a,M,'fs',fs);
%colormap(graycm);
ylim([0,10e3]);
shg

F = frame('dgtreal',g,a,M,'timeinv');
Fd = framedual(F);
maxit = 500;
cframe = frana(F,f);
[fhatGLA,relresglazero,~,chat] = frsynabs(F,abs(cframe),'rand','maxit',maxit,'fgriflim','tol',Eheapint,'Fd',Fd,'print');
nextprojc = frana(F,fhatGLA);
Egla = magnitudeerrdb(s,framecoef2native(F,nextprojc))

figure(2);clf;
plotdgtrealphasediff(angle(framecoef2native(F,chat)),angle(c),s,thrplot,a,M,'fs',fs);
%colormap(graycm);
ylim([0,10e3]);
shg

figure(10);clf;
logs = log(s);
maxlogs = max(logs(:));
logs(logs<maxlogs-11) = maxlogs - 11;
plotdgtreal(log(s),a,M,fs,'lin','clim',[-11,maxlogs]);
ylim([0,10000]);
%colormap(graycm);