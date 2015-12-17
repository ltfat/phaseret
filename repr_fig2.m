function repr_fig2

basepath = fileparts(which(mfilename));
exportdir = [basepath,filesep,'texexport'];

if ~exist(exportdir,'dir')
    mkdir(exportdir);
end

[f,fs] = gspi;
f = f(1:4*2048);

a = 16;
M = 2048;
db = 50;

g = 'gauss'; 

gamma = 0.9;
graycm = flipud(gray.^gamma);

thr = 10^(-db/20);

c = dgtreal(f,g,a,M,'timeinv');
s = abs(c);


[chat] = spsireal(abs(c),a,M,'timeinv');
fhatSpsi = idgtreal(chat,{'dual',g},a,M,numel(f),'timeinv');
Espsi = magnitudeerrdb(s,dgtreal(fhatSpsi,g,a,M,'timeinv'))

figure(2);
plotdgtrealphasediff(angle(c),angle(chat),abs(c),thr,a,M,'fs',fs);
colormap(graycm);
ylim([0,10e3]);
shg
saveas(gcf,[exportdir,filesep,'gspi_phasediffspsi.png']);


[chat,onset] = phaserecunwrapreal(c,g,a,M,'timeinv');
fhatMagron = idgtreal(chat,{'dual',g},a,M,numel(f),'timeinv');
Emagron = magnitudeerrdb(s,dgtreal(fhatMagron,g,a,M,'timeinv'))

figure(1);clf;
plotdgtrealphasediff(angle(c),angle(chat),abs(c),thr,a,M,'fs',fs);
colormap(graycm);
ylim([0,10e3]);
shg
saveas(gcf,[exportdir,filesep,'gspi_phasediffunwrap.png']);



chat = constructphasereal(abs(c),g,a,M,thr,'timeinv');
fhatHeapint = idgtreal(chat,{'dual',g},a,M,numel(f),'timeinv');
Eheapint = magnitudeerrdb(s,dgtreal(fhatHeapint,g,a,M,'timeinv'))

figure(3);clf;
plotdgtrealphasediff(angle(c),angle(chat),s,thr,a,M,'fs',fs);
colormap(graycm);
ylim([0,10e3]);
shg
saveas(gcf,[exportdir,filesep,'gspi_phasediffheapint.png']);

figure(4);clf;
plotdgtreal(s,a,M,fs,'dynrange',-20*log10(thr));
colormap(graycm);
ylim([0,10e3]);
saveas(gcf,sprintf([exportdir,filesep,'gspi_%i.png'],numel(f)));

