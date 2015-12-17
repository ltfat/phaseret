function comparetoall(databasePath,subdirs,varargin)


definput.keyvals.maxit=200;
definput.keyvals.maxsamples=[];
definput.keyvals.a=256;
definput.keyvals.M=1024;
definput.keyvals.thr1=1e-1;
definput.keyvals.thr2=1e-4;
definput.keyvals.exportdir = [];
definput.keyvals.expname = 'experiment';
definput.keyvals.rtisiiters = [1,2,3,4,5,12,25];
definput.flags.winmod = {'gauss','truncgauss','hann','hamming'};
definput.flags.writewavs = {'donotstorewavs','storewavs'};
definput.flags.warmstart = {'nowarmstart','wsheapint','wsspsi'};
definput.flags.lbfgs = {'dolbfgs','excludelbfgs'};
definput.keyvals.maxwavno = [];
[flags,kv,a,M]=ltfatarghelper({'a','M'},definput,varargin);
maxit = kv.maxit;


numFiles = 0;


prefix = [databasePath,filesep];
speakers = subdirs;


CheapintAll = 0;
CspsiAll = 0;
glarandAll = zeros(maxit,1);
glazeroAll = zeros(maxit,1);
fglaAll = zeros(maxit,1);
lbfgsAll = zeros(maxit,1);
leglaAll = zeros(maxit,1);
fleglaAll = zeros(maxit,1);
rtisiAll = zeros(numel(kv.rtisiiters),1);

lbfgsInvalidFilesNo = 0;

for jj = 1:numel(speakers)
    speaker = speakers{jj};
    
    dirname = [prefix,speaker,filesep];
    allwavs = dir([dirname,'*.wav']);
    
    if ~isempty(kv.maxwavno)
        allwavs = allwavs(1:min([kv.maxwavno,numel(allwavs)]));
    end
    
    numFiles = numFiles + numel(allwavs);
    
    allwavsCell = arrayfun(@(wav) [dirname,wav.name],allwavs,'UniformOutput',0);
    
    for ii=1:numel(allwavs)
        
        tic;
        
        Cres = compareonewav(allwavsCell{ii},kv,flags);
        
        toc;
        
        [Cheapint,Cspsi,Crelres,Cfinal] = deal(Cres{:});
        
        Cfinaldb = cellfun(@(cEl) 20*log10(cEl),Cfinal,'UniformOutput',0);
        Cheapintdb = 20*log10(Cheapint);
        Cspsidb = 20*log10(Cspsi);
        
            
        if ~flags.do_nowarmstart
            [relresglazero,relresfgla,relreslbfgs,relreslegla,relresflegla] = deal(Crelres{:});
            [Cglazero,Cfastgla,Clbfgs,Clegla,Cflegla] = deal(Cfinaldb{:});
            
            if ~flags.do_excludelbfgs
                fprintf(['Heapint: %.2f, SPSI: %.2f, GLA: %.2f,',...
                    'FGLA: %.2f, lBFGS: %.2f, leLGA: %.2f, FleLGA: %.2f, file: %s \n'],...
                    Cheapintdb,Cspsidb,Cglazero,Cfastgla,Clbfgs,Clegla,Cflegla,allwavsCell{ii});
            else
                fprintf(['Heapint: %.2f, SPSI: %.2f, GLA: %.2f,',...
                    'FGLA: %.2f, leLGA: %.2f, FleLGA: %.2f, file: %s \n'],...
                    Cheapintdb,Cspsidb,Cglazero,Cfastgla,Clegla,Cflegla,allwavsCell{ii});
            end
        else
            [relresglazero,relresglarand,relresfgla,relreslbfgs,relreslegla,relresflegla,relresrtisi]=deal(Crelres{:});
            [Cglazero,Cglarand,Cfastgla,Clbfgs,Clegla,Cflegla,Crtisi] = deal(Cfinaldb{:});
            
            if ~flags.do_excludelbfgs
                fprintf(['Heapint: %.2f, SPSI: %.2f, GLA rand: %.2f, GLA zero: %.2f,',...
                    'FGLA: %.2f, lBFGS: %.2f, leLGA: %.2f, FleLGA: %.2f, RTISI: %.2f, file: %s \n'],...
                    Cheapintdb,Cspsidb,Cglarand,Cglazero,Cfastgla,Clbfgs,Clegla,Cflegla,Crtisi,allwavsCell{ii});
            else
                fprintf(['Heapint: %.2f, SPSI: %.2f, GLA rand: %.2f, GLA zero: %.2f,',...
                    'FGLA: %.2f, leLGA: %.2f, FleLGA: %.2f, RTISI: %.2f, file: %s \n'],...
                    Cheapintdb,Cspsidb,Cglarand,Cglazero,Cfastgla,Clegla,Cflegla,Crtisi,allwavsCell{ii});                
            end
        end

        
        
        
        CheapintAll = CheapintAll + Cheapint;
        CspsiAll = CspsiAll + Cspsi;
        glazeroAll = glazeroAll + postpad(relresglazero,maxit,relresglazero(end));
        fglaAll = fglaAll + postpad(relresfgla,maxit,relresfgla(end));
        if ~isnan(relreslbfgs)
            lbfgsAll = lbfgsAll + postpad(relreslbfgs,maxit,relreslbfgs(end));
        else
            lbfgsInvalidFilesNo = lbfgsInvalidFilesNo + 1;
        end
        leglaAll = leglaAll + postpad(relreslegla,maxit,relreslegla(end));
        fleglaAll = fleglaAll + postpad(relresflegla,maxit,relresflegla(end));
        
        if flags.do_nowarmstart
            glarandAll = glarandAll + postpad(relresglarand,maxit,relresglarand(end));
        end
        
        
        rtisiAll = rtisiAll + relresrtisi;
        
        
        
    end
end


figure(10);

glazeroMean = glazeroAll./numFiles;
fglaMean = fglaAll./numFiles;
lbfgsMean = lbfgsAll./(numFiles-lbfgsInvalidFilesNo);
leglaMean = leglaAll./numFiles;
fleglaMean = fleglaAll./numFiles;
rtisiMean = rtisiAll./numFiles;
rtisiXaxis = kv.rtisiiters*(M/a);
CspsiAllMean = CspsiAll./numFiles;
CheapintAllMean = CheapintAll./numFiles;

if ~flags.do_nowarmstart
    if ~flags.do_excludelbfgs
        plot(20*log10([fglaMean,lbfgsMean,leglaMean,fleglaMean,glazeroMean]));
        legend({'FGLA','lBFGS','leGLA','FleGLA','GLA'})
    else
        plot(20*log10([fglaMean,leglaMean,fleglaMean,glazeroMean]));
        legend({'FGLA','leGLA','FleGLA','GLA'})        
    end
else
    glarandMean = glarandAll./numFiles;
    if ~flags.do_excludelbfgs
        plot(20*log10([fglaMean,lbfgsMean,leglaMean,fleglaMean,glazeroMean,glarandMean]));
        hold on; plot(rtisiXaxis,20*log10(rtisiMean),'ko--');hold off;
        legend({'FGLA','lBFGS','leGLA','FleGLA','GLA zero','GLA rand','RTISI'})
    else
        plot(20*log10([fglaMean,leglaMean,fleglaMean,glazeroMean,glarandMean]));
        hold on; plot(rtisiXaxis,20*log10(rtisiMean),'ko--');hold off;
        legend({'FGLA','leGLA','FleGLA','GLA zero','GLA rand','RTISI'})
    end
end

if flags.do_wsspsi
    h = line([1,maxit],20*log10([CspsiAllMean,CspsiAllMean]));
else
    h = line([1,maxit],20*log10([CheapintAllMean,CheapintAllMean]));
end
set(h,'LineStyle','--');
set(h,'Color','k');

if ~isempty(kv.exportdir)
    wsstring = '';
    if ~flags.do_nowarmstart
        wsstring = flags.warmstart;
        Cres = {CheapintAllMean,CspsiAllMean,fglaMean,lbfgsMean,leglaMean,fleglaMean,glazeroMean};
        save(sprintf([kv.exportdir,filesep,'%scomparison_%s_%s.mat'],flags.winmod,kv.expname,wsstring),...
        'Cres');
    else
        Cres = {CheapintAllMean,CspsiAllMean,fglaMean,lbfgsMean,leglaMean,fleglaMean,glazeroMean,glarandMean};
        save(sprintf([kv.exportdir,filesep,'%scomparison_%s_%s.mat'],flags.winmod,kv.expname,wsstring),...
        'Cres');
    end
    

    saveas(gcf,sprintf([kv.exportdir,filesep,'%scomparison_%s_%s.fig'],flags.winmod,kv.expname,wsstring),'fig');
    saveas(gcf,sprintf([kv.exportdir,filesep,'%scomparison_%s_%s.png'],flags.winmod,kv.expname,wsstring),'png');
end


function [Cres] = compareonewav(wavfile,kv,flags)
wsstring = '';
if ~flags.do_nowarmstart
    wsstring = flags.warmstart;
end

a = kv.a;
M = kv.M;
thr1 = kv.thr1;
thr2 = kv.thr2;
maxit = kv.maxit;
[~,name] = fileparts(wavfile);

[f,fs] = wavreadwrapper(wavfile);
f = normalize(f(:,1),'wav');

if ~isempty(kv.maxsamples)
    f = f(1:min([kv.maxsamples,numel(f)]));
end

Ls= numel(f);

if flags.do_storewavs
    outputdir = [kv.exportdir,filesep,mfilename,'_',kv.expname,'_',flags.winmod];
    nameoutdir = [outputdir,filesep,name,filesep];
    
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    end
    
    if ~exist(nameoutdir, 'dir')
        mkdir(nameoutdir);
    end
    
    wavsave(normalize(f,'wav'),fs,[nameoutdir,name,'.wav']);
end

L = dgtlength(Ls,a,M);
tfr = a*M/L;
ggauss = {'gauss',tfr,'inf'};
f = postpad(f,L);

switch flags.winmod
    case 'gauss'
        g = ggauss;
        gd = gabdual(g,a,M,L);
    case 'truncgauss'
        gnum = gabwin(ggauss,a,M,L);
        g = long2fir(gnum,M);
        gd = gabdual(g,a,M,L);
    case 'hann'
        gl = M;
        atheight = 0.277;
        g = {'hann',gl};
        gnum = gabwin(g,a,M);
        w = winwidthatheight(gnum,atheight);
        ggauss = {'gauss','width',w,'atheight',atheight,'inf'};
        gd = gabdual(g,a,M,L);
    case 'hamming'
        gl = M;
        atheight = 0.306;
        g = {'hamming',gl};
        gnum = gabwin(g,a,M);
        w = winwidthatheight(gnum,atheight);
        ggauss = {'gauss','width',w,'atheight',atheight,'inf'};
        gd = gabdual(g,a,M,L);
end

c = dgtreal(f,g,a,M);

s = abs(c);
F = frame('dgtreal',g,a,M);
Fd = frame('dgtreal',gd,a,M);
cframe = framenative2coef(F,c);

[chatint,newphase,usedmask] = constructphasereal(abs(c),ggauss,a,M,thr1);
chatint = constructphasereal(chatint,ggauss,a,M,thr2,'mask',usedmask,'usephase',newphase);

chatintframe = framenative2coef(F,chatint);
fhat = idgtreal(chatint,gd,a,M);

%%
%
disp('Proposed:');
nextprojc = dgtreal(fhat,g,a,M);
Cheapint = magnitudeerr(c,nextprojc);

if flags.do_storewavs
    os = [nameoutdir,name,'_heapint'];
    writetex([os,'.tex'],20*log10(Cheapint));
    wavsave(normalize(fhat,'wav'),fs,[os,'.wav']);
end

disp('SPSI:');
chatspsi = spsireal(abs(c),a,M);
chatspsiframe = framenative2coef(F,chatspsi);
fhatspsi = idgtreal(chatspsi,gd,a,M);
nextprojc = dgtreal(fhatspsi,g,a,M);
Cspsi = magnitudeerr(c,nextprojc);

if flags.do_storewavs
    os = [nameoutdir,name,'_spsi'];
    writetex([os,'.tex'],20*log10(Cspsi));
    wavsave(normalize(fhatspsi,'wav'),fs,[os,'.wav']);
end

if flags.do_wsheapint
    c = chatint;
    cframe = chatintframe;
elseif flags.do_wsspsi
    c = chatspsi;
    cframe = chatspsiframe;
else
    c = abs(c);
    cframe = abs(cframe);
end

%%
%
if flags.do_nowarmstart
    disp('Grif-Lim rand:');
    [fhat,relresglarand,~,~] = frsynabs(F,cframe,'maxit',maxit,'griflim','tol',1e-10,'quiet','rand','Fd',Fd);
    nextprojc = frana(F,fhat);
    Cglarand = magnitudeerr(c,framecoef2native(F,nextprojc));
    
    if flags.do_storewavs
        os = [nameoutdir,name,'_GLArand',wsstring];
        writetex([os,'.tex'],20*log10(Cglarand));
        wavsave(normalize(fhat,'wav'),fs,[os,'.wav']);
    end
end

disp('Grif-Lim:');
[fhat,relresglazero,~,chat] = frsynabs(F,cframe,'maxit',maxit,'griflim','tol',1e-10,'quiet','Fd',Fd);

nextprojc = frana(F,fhat);
Cglazero = magnitudeerr(c,framecoef2native(F,nextprojc));

if flags.do_storewavs
    os = [nameoutdir,name,'_GLArand',wsstring];
    writetex([os,'.tex'],20*log10(Cglazero));
    wavsave(normalize(fhat,'wav'),fs,[os,'.wav']);
end
%%
%
disp('Fast Grif-Lim')
[fhat,relresfgla,~,chat] = frsynabs(F,cframe,'maxit',maxit,'fgriflim','tol',1e-10,'quiet','Fd',Fd);

nextprojc = frana(F,fhat);
Cfastgla = magnitudeerr(c,framecoef2native(F,nextprojc));

if flags.do_storewavs
    os = [nameoutdir,name,'_FGLA',wsstring];
    writetex([os,'.tex'],20*log10(Cfastgla));
    wavsave(normalize(fhat,'wav'),fs,[os,'.wav']);
end


if ~flags.do_excludelbfgs
    disp('lBFGS');
    % minFunc migth not finish properly.
    try
        pause off;
        [fhat,relreslbfgs] = frsynabs(F,cframe,'maxit',maxit,'bfgs','p',2/3,'tol',1e-9,'quiet','Fd',Fd);
        
        nextprojc = frana(F,fhat);
        Clbfgs = magnitudeerr(c,framecoef2native(F,nextprojc));
        
        if flags.do_storewavs
            os = [nameoutdir,name,'_lBFGS',wsstring];
            writetex([os,'.tex'],20*log10(Clbfgs));
            wavsave(normalize(fhat,'wav'),fs,[os,'.wav']);
        end
    catch
        err = lasterror;
        tocompare = 'frsynabs';
        if strcmpi(err.stack(1).name,tocompare)
            fprintf(err.message);
            % minFunc failed for some reason, we do not want this result to be
            % included in the average.
            relreslbfgs = nan;
            Clbfgs = nan;
            
            if flags.do_storewavs
                os = [nameoutdir,name,'_lBFGS',wsstring];
                writetex([os,'.tex'],Clbfgs);
            end
            
        else
            pause on;
            rethrow(err);
        end
    end
    pause on;
else
    Clbfgs = nan;
    relreslbfgs = nan;
end



disp('leGLA');
[fhat,relreslegla] = leglareal(c,g,a,M,'maxit',maxit,'tol',1e-9,'modtrunc','onthefly');

nextprojc = frana(F,fhat);
Clegla = magnitudeerr(c,framecoef2native(F,nextprojc));

if flags.do_storewavs
    os = [nameoutdir,name,'_leGLA',wsstring];
    writetex([os,'.tex'],20*log10(Clegla));
    wavsave(normalize(fhat,'wav'),fs,[os,'.wav']);
end

disp('FleGLA');
[fhat,relresflegla] = leglareal(c,g,a,M,'maxit',maxit,'tol',1e-9,'modtrunc','onthefly','flegla');

nextprojc = frana(F,fhat);
Cflegla = magnitudeerr(c,framecoef2native(F,nextprojc));

if flags.do_storewavs
    os = [nameoutdir,name,'_FleGLA',wsstring];
    writetex([os,'.tex'],20*log10(Cflegla));
    wavsave(normalize(fhat,'wav'),fs,[os,'.wav']);
end

if flags.do_nowarmstart
    relresrtisi = zeros(numel(kv.rtisiiters),1);
    
    disp('RTISI');
    for iIt=1:numel(kv.rtisiiters)
        [fhat,relresrtisi(iIt)] = rtisireal(c,g,a,M,'maxit',kv.rtisiiters(iIt));
    end

    nextprojc = frana(F,fhat);
    Crtisi = magnitudeerr(c,framecoef2native(F,nextprojc));

    if flags.do_storewavs
        os = [nameoutdir,name,'_RTISILA',wsstring];
        writetex([os,'.tex'],20*log10(Crtisi));
        wavsave(normalize(fhat,'wav'),fs,[os,'.wav']);
    end
end


if ~flags.do_nowarmstart
    Crelres = {relresglazero,relresfgla,relreslbfgs,relreslegla,relresflegla};
    Cfinal = {Cglazero,Cfastgla,Clbfgs,Clegla,Cflegla};
else
    Crelres = {relresglazero,relresglarand,relresfgla,relreslbfgs,relreslegla,relresflegla,relresrtisi};
    Cfinal = {Cglazero,Cglarand,Cfastgla,Clbfgs,Clegla,Cflegla,Crtisi};
end

Cres = {Cheapint,Cspsi,Crelres,Cfinal};

function writetex(s,val)
fileID = fopen(s,'w');
fprintf(fileID,'%.2f\n',val);
fclose(fileID);


