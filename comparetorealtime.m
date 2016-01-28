function comparetorealtime(databasePath,subdirs,varargin)

definput.keyvals.a=256;
definput.keyvals.M=2048;
definput.keyvals.maxsamples=[];
definput.keyvals.gl=[];
definput.keyvals.thr1=1e-1;
definput.keyvals.thr2=1e-10;
definput.keyvals.perframeit=8;
definput.keyvals.exportdir = [];
definput.keyvals.expname = 'experiment';
definput.flags.winmod = {'truncgauss','hann','hamming','blackman'};
definput.flags.writewavs = {'donotstorewavs','storewavs'};
[flags,kv,a,M]=ltfatarghelper({'a','M'},definput,varargin);
thr1 = kv.thr1;
thr2 = kv.thr2;

if isempty(kv.gl)
    kv.gl = kv.M;
end

prefix = [databasePath,filesep];

CheapintAll = [];
CspsiAll = [];
CrtheapintAll = [];
Crtheapint0All = [];
Crtisila0All = [];
Crtisila1All = [];

for jj= 1:numel(subdirs)
    subdir = subdirs{jj};
    
    dirname = [prefix,subdir,filesep];
    allwavs = dir([dirname,'*.wav']);
    
    for ii=1:numel(allwavs)
        wavfile = allwavs(ii).name;
        name = wavfile(1:strfind(wavfile,'.')-1);
        outputdir = [kv.exportdir,filesep,mfilename,'_',kv.expname,'_',flags.winmod];
        nameoutdir = [outputdir,filesep,name,filesep];
        
        if ~exist(outputdir, 'dir')
            mkdir(outputdir);
        end
        
        if ~exist(nameoutdir, 'dir')
            mkdir(nameoutdir);
        end
        
        [f,fs] = wavreadwrapper([dirname,wavfile]);
        f = normalize(f(:,1),'wav');
        if ~isempty(kv.maxsamples)
            f = f(1:min([kv.maxsamples,numel(f)]));
        end
        Ls = numel(f);
        
        if flags.do_storewavs
            wavsave(normalize(f,'wav'),fs,[nameoutdir,name,'.wav']);
        end
        
        L = dgtlength(Ls,a,M);
        f = postpad(f,L);
        
        switch flags.winmod
            case 'truncgauss'
                h = 0.0001;
                gnum = gabwin({'gauss','width',kv.gl,'atheight',h},a,M,10*M);
                gnum = long2fir(gnum,kv.gl);
                lambda = -pi/4*(kv.gl)^2/log(h);
            case 'hann'
                g = {'hann',kv.gl};
                gnum = gabwin(g,a,M);
                lambda = 0.2562*(kv.gl)^2;
            case 'hamming'
                g = {'hamming',kv.gl};
                gnum = gabwin(g,a,M);
                lambda = 0.29765*(kv.gl)^2;
            case 'blackman'
                g = {'blackman',kv.gl};
                gnum = gabwin(g,a,M);
                lambda = 0.17937*(kv.gl)^2;               
        end
        gdnum = gabdual(gnum,a,M,L);
        
        
        c = comp_sepdgtreal(f,gnum,a,M,1);
        
        s = abs(c);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        chatint = constructphasereal(s,{'gauss',lambda/L},a,M,[thr1,thr2],'timeinv');
        HEAPINTtime = toc;
        fhatint = comp_isepdgtreal(chatint,gdnum,L,a,M,1);
        
        nextprojc = dgtreal(fhatint,gnum,a,M,'timeinv');
        Cheapint = magnitudeerr(c,nextprojc);
        Cheapintdb = 20*log10(Cheapint);
        CheapintAll(end+1) = Cheapint;
        
        if flags.do_storewavs
            writetex(sprintf([nameoutdir,name,'_heapint%s%s.tex'],flags.winmod,kv.expname),Cheapintdb);
            wavsave(normalize(fhatint,'wav'),fs,[nameoutdir,name,'_heapint.wav']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        tic;
        chatspsi = spsireal(s,a,M,'timeinv');
        SPSItime = toc;
        fhatspsi = idgtreal(chatspsi,gdnum,a,M,'timeinv');
        
        nextprojc = dgtreal(fhatspsi,gnum,a,M,'timeinv');
        Cspsi = magnitudeerr(c,nextprojc);
        Cspsidb = 20*log10(Cspsi);
        CspsiAll(end+1) = Cspsi;
        
        if flags.do_storewavs
            writetex(sprintf([nameoutdir,name,'_spsi%s%s.tex'],flags.winmod,kv.expname),Cspsidb);
            wavsave(normalize(fhatspsi,'wav'),fs,[nameoutdir,name,'_spsi.wav']);
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        chatrtint = simrealtime(s,lambda,a,M,[thr1,thr2],'timeinv');
        HEAPINTRTtime = toc;
        fhatrtint = comp_isepdgtreal(chatrtint,gdnum,L,a,M,1);
        
        nextprojc = dgtreal(fhatrtint,gnum,a,M,'timeinv');
        Crtheapint = magnitudeerr(c,nextprojc);
        Crtheapintdb = 20*log10(Crtheapint);
        CrtheapintAll(end+1) = Crtheapint; 
        
        if flags.do_storewavs
            writetex(sprintf([nameoutdir,name,'_rtheapint%s%s.tex'],flags.winmod,kv.expname),Crtheapintdb);
            wavsave(normalize(fhatrtint,'wav'),fs,[nameoutdir,name,'_rtheapint.wav']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        chatrtint0 = simrealtime(s,lambda,a,M,[thr1,thr2],'timeinv','causal');
        HEAPINTRT0time = toc;
        fhatrtint0 = comp_isepdgtreal(chatrtint0,gdnum,L,a,M,1);
        
        nextprojc = dgtreal(fhatrtint0,gnum,a,M,'timeinv');
        Crtheapint0 = magnitudeerr(c,nextprojc);
        Crtheapint0db = 20*log10(Crtheapint0);
        Crtheapint0All(end+1) = Crtheapint0;
         
        if flags.do_storewavs
            writetex(sprintf([nameoutdir,name,'_rtheapint0%s%s.tex'],flags.winmod,kv.expname),Crtheapint0db);
            wavsave(normalize(fhatrtint0,'wav'),fs,[nameoutdir,name,'_rtheapint0.wav']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        chatrtisila0 = rtisirealforreal(s,gnum,a,M,'maxit',kv.perframeit,'lookahead',0,'timeinv');
        RTISILA0time = toc;
        fhatrtisila0 = comp_isepdgtreal(chatrtisila0,gdnum,L,a,M,1);
        
        nextprojc = dgtreal(fhatrtisila0,gnum,a,M,'timeinv');
        Crtisila0 = magnitudeerr(c,nextprojc);
        Crtisila0db = 20*log10(Crtisila0);
        Crtisila0All(end+1) = Crtisila0;
        
        if flags.do_storewavs
            writetex(sprintf([nameoutdir,name,'_rtisila0%s%s.tex'],flags.winmod,kv.expname),Crtisila0db);
            wavsave(normalize(fhatrtisila0,'wav'),fs,[nameoutdir,name,'_rtisila0.wav']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        chatrtisila1 = rtisirealforreal(s,gnum,a,M,'maxit',floor(kv.perframeit/2),'lookahead',1,'timeinv');
        RTISILA1time = toc;
        fhatrtisila1 = comp_isepdgtreal(chatrtisila1,gdnum,L,a,M,1);
        
        nextprojc = dgtreal(fhatrtisila1,gnum,a,M,'timeinv');
        Crtisila1 = magnitudeerr(c,nextprojc);
        Crtisila1db = 20*log10(Crtisila1);
        Crtisila1All(end+1) = Crtisila1;
        
        if flags.do_storewavs
            writetex(sprintf([nameoutdir,name,'_rtisila1%s%s.tex'],flags.winmod,kv.expname),Crtisila1db);
            wavsave(normalize(fhatrtisila1,'wav'),fs,[nameoutdir,name,'_rtisila1.wav']);
        end
        
        fprintf(['Heapint: %.2f dB, SPSI: %.2f dB, ',...
                 'RTheapint1: %.2f dB, RTISI1 %.2f dB, ',...
                 'RTheapint0: %.2f dB, RTISI0 %.2f dB',...
                 ' %s, window: %s \n'],...
                 Cheapintdb,Cspsidb,...
                 Crtheapintdb,Crtisila1db,...
                 Crtheapint0db,Crtisila0db,...
                 wavfile,flags.winmod);
        
    end
    
end


CheapintAllMeanDB = 20*log10(mean(CheapintAll));
CspsiAllMeanDB = 20*log10(mean(CspsiAll));
CrtheapintAllMeanDB = 20*log10(mean(CrtheapintAll));
Crtheapint0AllMeanDB = 20*log10(mean(Crtheapint0All));
Crtisila0AllMeanDB = 20*log10(mean(Crtisila0All));
Crtisila1AllMeanDB = 20*log10(mean(Crtisila1All));

disp('--------------------------MEAN------------------------------');
        fprintf(['Heapint: %.2f dB, SPSI: %.2f dB, ',...
                 'RTheapint1: %.2f dB, RTISI1 %.2f dB, ',...
                 'RTheapint0: %.2f dB, RTISI0 %.2f dB',...
                 ' window: %s \n'],...
                 CheapintAllMeanDB,CspsiAllMeanDB,...
                 CrtheapintAllMeanDB,Crtisila1AllMeanDB,...
                 Crtheapint0AllMeanDB,Crtisila0AllMeanDB,...
                 flags.winmod);
disp('------------------------------------------------------------');

if flags.do_storewavs
    writetex(sprintf([nameoutdir,'avg_heapint%s%s.tex'],flags.winmod,kv.expname),CheapintAllMeanDB);
    writetex(sprintf([nameoutdir,'avg_spsi%s%s.tex'],flags.winmod,kv.expname),CspsiAllMeanDB);
    writetex(sprintf([nameoutdir,'avg_rtheapint%s%s.tex'],flags.winmod,kv.expname),CrtheapintAllMeanDB);
    writetex(sprintf([nameoutdir,'avg_rtheapint0%s%s.tex'],flags.winmod,kv.expname),Crtheapint0AllMeanDB);
    writetex(sprintf([nameoutdir,'avg_rtisila0%s%s.tex'],flags.winmod,kv.expname),Crtisila0AllMeanDB);
    writetex(sprintf([nameoutdir,'avg_rtisila1%s%s.tex'],flags.winmod,kv.expname),Crtisila1AllMeanDB);
end


function writetex(s,val)
fileID = fopen(s,'w');
fprintf(fileID,'%.2f\n',val);
fclose(fileID);







