function comparetospsi(databasePath,subdirs,varargin)


definput.keyvals.a=256;
definput.keyvals.M=1024;
definput.keyvals.maxsamples=[];
definput.keyvals.gl=[];
definput.keyvals.thr1=1e-1;
definput.keyvals.thr2=1e-10;
definput.keyvals.exportdir = [];
definput.keyvals.expname = 'experiment';
definput.flags.winmod = {'gauss','truncgauss','hann','hamming'};
definput.flags.writewavs = {'donotstorewavs','storewavs'};
[flags,kv,a,M]=ltfatarghelper({'a','M'},definput,varargin);
thr1 = kv.thr1;
thr2 = kv.thr2;

if isempty(kv.gl)
    kv.gl = kv.M;
end

basepath = fileparts(which(mfilename));

prefix = [databasePath,filesep];
speakers = subdirs;

CheapintAll = [];
CspsiAll = [];
CunwrapAll = [];

for jj = 1:numel(speakers)
    speaker = speakers{jj};
    
    dirname = [prefix,speaker,filesep];
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
        
        tfr = a*M/L;
        ggauss = {'gauss',tfr,'inf'};
        
        switch flags.winmod
            case 'gauss'
                g = ggauss;
                gnum = gabwin(g,a,M,L);
            case 'truncgauss'
                gnum = gabwin(ggauss,a,M,L);
                g = long2fir(gnum,kv.gl);
                gnum = g;
            case 'hann'
                gl = kv.gl;
                atheight = 0.277;
                g = {'hann',gl};
                gnum = gabwin(g,a,M);
                w = winwidthatheight(gnum,atheight);
                ggauss = {'gauss','width',w,'atheight',atheight,'inf'};
            case 'hamming'
                gl = kv.gl;
                atheight = 0.306;
                g = {'hamming',gl};
                gnum = gabwin(g,a,M);
                w = winwidthatheight(gnum,atheight);
                ggauss = {'gauss','width',w,'atheight',atheight,'inf'};
        end
        gdnum = gabdual(gnum,a,M,L);
        
        % DO NOT CALL THIS FUNCTION DIRECTLY!! USE
        % c = dgtreal(f,gnum,a,M,'timeinv');
        % INSTEAD!!!!
         tic;
         c = comp_sepdgtreal(f,gnum,a,M,1);
         DGTtime = toc;
        
        s = abs(c);
        tic;
        cupd = s.*exp(1i*angle(c));
        UPDtime = toc;
        
        lala = cupd;
        
        tic;
        chatint = constructphasereal(s,ggauss,a,M,[thr1,thr2],'timeinv');
        HEAPINTtime = toc;
        
        
        %fhat = idgtreal(chatint,gdnum,a,M,'timeinv');
        tic;
        fhat = comp_isepdgtreal(chatint,gdnum,L,a,M,1);
        IDGTtime = toc;
        
        if flags.do_storewavs
            wavsave(normalize(fhat,'wav'),fs,[nameoutdir,name,'_heapint.wav']);
        end
        
        tic;
        chatspsi = spsireal(s,a,M,'timeinv');
        SPSItime = toc;
        fhatspsi = idgtreal(chatspsi,gdnum,a,M,'timeinv');
        
        if flags.do_storewavs
            wavsave(normalize(fhatspsi,'wav'),fs,[nameoutdir,name,'_spsi.wav']);
        end
        
%         chatunwrap = phaserecunwrapreal(s,g,a,M);
%         fhatunwrap = idgtreal(chatunwrap,{'dual',g},a,M);
%         
%         if flags.do_storewavs
%             wavsave(normalize(fhatunwrap,'wav'),fs,[nameoutdir,name,'_unwrap.wav']);
%         end       
        
        
        nextprojc = dgtreal(fhat,gnum,a,M,'timeinv');
        Cheapint = magnitudeerr(c,nextprojc);
        Cheapintdb = 20*log10(Cheapint);
        
        nextprojc = dgtreal(fhatspsi,gnum,a,M,'timeinv');
        Cspsi = magnitudeerr(c,nextprojc);
        Cspsidb = 20*log10(Cspsi);
        
%         nextprojc = dgtreal(fhatunwrap,g,a,M);
%         Cunwrap = magnitudeerr(c,nextprojc);
%         Cunwrapdb = 20*log10(Cunwrap);
        
        if flags.do_storewavs
            fileID = fopen(sprintf([nameoutdir,name,'_heapint%s%s.tex'],flags.winmod,kv.expname),'w');
            fprintf(fileID,'%.2f\n',Cheapintdb);
            fclose(fileID);
            fileID = fopen(sprintf([nameoutdir,name,'_spsi%s%s.tex'],flags.winmod,kv.expname),'w');
            fprintf(fileID,'%.2f\n',Cspsidb);
            fclose(fileID);
%             fileID = fopen(sprintf([nameoutdir,name,'_unwrap%s%s.tex'],flags.winmod,kv.expname),'w');
%             fprintf(fileID,'%.2f\n',Cunwrapdb);
%             fclose(fileID);
        end
        
        CheapintAll(end+1) = Cheapint;
        CspsiAll(end+1) = Cspsi;
%        CunwrapAll(end+1) = Cunwrap;
        
        %fprintf('Heapint: %.2f, SPSI: %.2f, UNWRAP: %.2f, %s , window: %s \n',Cheapintdb,Cspsidb,Cunwrapdb,wavfile,flags.winmod);
        fprintf(['Heapint: %.2f dB, time: %.4fs, SPSI: %.2f dB, ',...
                 'time %.4fs, it. time %.4fs, %s, window: %s \n'],...
                 Cheapintdb,HEAPINTtime,Cspsidb,SPSItime,...
                 DGTtime+IDGTtime + UPDtime,wavfile,flags.winmod);
      
        
        
    end
end

CheapintAllMeanDB = 20*log10(mean(CheapintAll));
CspsiAllMeanDB = 20*log10(mean(CspsiAll));
%CunwrapAllMeanDB = 20*log10(mean(CunwrapAll));

%fprintf('Mean Heapint: %.2f, Mean SPSI: %.2f,  Mean UNWRAP: %.2f\n',CheapintAllMeanDB,CspsiAllMeanDB,CunwrapAllMeanDB);
fprintf('Mean Heapint: %.2f, Mean SPSI: %.2f \n',CheapintAllMeanDB,CspsiAllMeanDB);


if flags.do_storewavs
    fileID = fopen(sprintf([nameoutdir,'avg_heapint%s%s.tex'],flags.winmod,kv.expname),'w');
    fprintf(fileID,'%.2f\n',CheapintAllMeanDB);
    fclose(fileID);
    fileID = fopen(sprintf([nameoutdir,'avg_spsi%s%s.tex'],flags.winmod,kv.expname),'w');
    fprintf(fileID,'%.2f\n',CspsiAllMeanDB);
    fclose(fileID);
%     fileID = fopen(sprintf([nameoutdir,'avg_unwrap%s%s.tex'],flags.winmod,kv.expname),'w');
%     fprintf(fileID,'%.2f\n',CunwrapAllMeanDB);
%     fclose(fileID);
end

if ~isempty(kv.exportdir)
    fileID = fopen(sprintf([kv.exportdir,filesep,'heapint%s%s.tex'],flags.winmod,kv.expname),'w');
    fprintf(fileID,'%.2f\n',CheapintAllMeanDB);
    fclose(fileID);
    fileID = fopen(sprintf([kv.exportdir,filesep,'spsi%s%s.tex'],flags.winmod,kv.expname),'w');
    fprintf(fileID,'%.2f\n',CspsiAllMeanDB);
    fclose(fileID);
%     fileID = fopen(sprintf([kv.exportdir,filesep,'unwrap%s%s.tex'],flags.winmod,kv.expname),'w');
%     fprintf(fileID,'%.2f\n',CunwrapAllMeanDB);
%     fclose(fileID);
end





