function repr_fig1
% This experiments show performance of the phase integration depending on
% the step size in time

% 
basepath = fileparts(which(mfilename));
exportdir = [basepath,filesep,'texexport'];

if ~exist(exportdir,'dir')
    mkdir(exportdir);
end


[f,fs] = greasy; Ls= numel(f); f = postpad(f,ceil(Ls/256)*256);
Ls= numel(f);

Marr = {Ls, Ls, Ls };
aarr = {1, 16 , 32 };
db = 60; 

for ii=1:numel(Marr)
    M = Marr{ii}; a = aarr{ii};
    L = dgtlength(Ls,a,M);

    g = {'gauss',1,'inf'};
    gnum = gabwin(g,a,M,L);

    [c,~,gnum] = dgtreal(f,gnum,a,M);
    s = abs(c);
    oldphase = angle(c);
    tol = 10^(-db/20);

    chatint = constructphasereal(abs(c),g,a,M,tol);
    [fhat,gd] = idgtreal(chatint,{'dual',gnum},a,M,Ls);

    nextprojc = dgtreal(fhat,gnum,a,M);

    C = magnitudeerrdb(s,nextprojc)

    fileID = fopen(sprintf([exportdir,filesep,'phasediffC_%i.tex'],ii),'w');
    fprintf(fileID,'%.2f\n',C);
    fclose(fileID);

    fileID = fopen(sprintf([exportdir,filesep,'phasediffa_%i.tex'],ii),'w');
    fprintf(fileID,'%d\n',a);
    fclose(fileID);

    cm = flipud(gray);

    figure(1);
    plotdgtrealphasediff(angle(chatint),oldphase,s,tol,a,M,'fs',fs);
    colormap(cm);

    saveas(gcf,sprintf([exportdir,filesep,'phasediff_%i.png'],ii));

    figure(2);
    plotdgtreal(c,a,M,'dynrange',db,'fs',fs);
    colormap(cm);

    saveas(gcf,sprintf([exportdir,filesep,'greasy_%i.png'],ii));
    
end

