[f,fs] = greasy;
Ls = numel(f);

a = 16;
M = 1024;
M2 = floor(M/2) + 1;
thr = 1e-1;
thr2 = 1e-10;%1e-10;

L=dgtlength(Ls,a,M);
b = L/M;
N = L/a;

ggauss = 'gauss';


c = dgtreal(f,ggauss,a,M,'timeinv');
cfreq = dgtreal(f,ggauss,a,M,'freqinv');
s = abs(c);

chat = constructphasereal(abs(c),ggauss,a,M,[thr,thr2],'timeinv');
HEAPINTtime = toc;

figure(1);
plotdgtreal(s,a,M,'fs',fs,'dynrange',80);
saveas(gcf,'images/title_greasy.png');

fhatHeapint = idgtreal(chat,{'dual',g},a,M,numel(f),'timeinv');
nextprojc = dgtreal(fhatHeapint,g,a,M,'timeinv');
EheapintF = magnitudeerrdb(s,nextprojc)

figure(2);
plotdgtrealphasediff(angle(c),angle(chat),s,thr2,a,M,'fs',fs);
%colormap(graycm);
ylim([0,8e3]);
shg

saveas(gcf,'images/title_greasy_diff.png');





