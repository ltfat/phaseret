function demo_blockproc_phaseret2(source,varargin) %RUNASSCRIPT
%DEMO_BLOCKPROC_PHASERET2 STFT phase reconstruction algorithms comparison
%   Usage: demo_blockproc_phaseret2('gspi.wav')
%          demo_blockproc_phaseret2(...,'a',a)` 
%
%   For additional help call |demo_blockproc_phaseret2| without arguments.
%
%   This file does real-time reconstruction from the STFT magnitude using
%   4 algorithms. No effect, just plain reconstruction.
%
%   The order of the algorithms (Algorithm slider) is the following::
%       0 - RTPGHI
%       1 - RTISI-LA
%       2 - GSRTISI-LA
%       3 - RTPGHI + GSRTISI-LA
%
%   Please note that you might experience some artifacts when switching
%   between the algorithms or when changing their parameters. This
%   does not have anything to do with the algorithms themselves.
%       
%   `demo_blockproc_phaseret(...,'a',a)` allows setting the hop factor.
%   Default value is 256.
%
%   .. image:: ../images/phaseret2_common.png
%   .. image:: ../images/phaseret2_rtpghi.png
%   .. image:: ../images/phaseret2_rtisila.png
%   .. image:: ../images/phaseret2_gsrtisila.png
%   .. image:: ../images/phaseret2_proposed.png
%
%   References: ltfatnote048 ltfatnote043
%

% AUTHOR: Zdenek Prusa

if demo_blockproc_header(mfilename,nargin)
    return;
end

[a,varargin] = parsevararginfor('a',varargin,256);
[M,varargin] = parsevararginfor('M',varargin,2048);
[gl,varargin] = parsevararginfor('gl',varargin,2048);

% Window support is M
h = 0.01;
g = gabwin({'gauss','width',gl,'atheight',0.01,'inf'},a,gl,10*gl);
g = gabwin(long2fir(g,gl),a,M);
gd = gabdual(g,a,M);

lookaheadmax = 11;

% Control pannel (Java object)
% Each entry determines one parameter to be changed during the main loop
% execution.
commonpar = {
    {'GdB','Gain',-20,20,0,21},...
    {'Bypass' ,'Bypass' ,0,1,0,2},...
    {'algno' ,'Algorithm' ,0,3,0,4},...
    };
... % Alg 0 params:
    alg0par = {
    {'Thr','Thr',0,16,6,100},...
    {'rtpghila0' ,'RTPGHI(N)' ,0,1,1,2},...
    };
alg0parNames = cellfun(@(aEl) aEl{1},alg0par,'UniformOutput',0);

alg1par = {
    {'lookahead1','Lookahead',0,lookaheadmax,2,lookaheadmax+1},...
    {'maxit1','Max. per-frame iter.',1,100,12,100},...
    };
alg1parNames = cellfun(@(aEl) aEl{1},alg1par,'UniformOutput',0);

alg2par = {
    {'lookahead2','Lookahead',0,lookaheadmax,2,lookaheadmax+1},...
    {'maxit2','Max. per-frame iter.',1,100,12,100},...
    };
alg2parNames = cellfun(@(aEl) aEl{1},alg2par,'UniformOutput',0);

alg3par = {
    {'lookahead3','Lookahead',0,lookaheadmax,2,lookaheadmax+1},...
    {'maxit3','Max. per-frame iter.',1,100,12,100},...
    {'rtpghila3' ,'RTPGHI(N)' ,0,1,1,2},...
    };
alg3parNames = cellfun(@(aEl) aEl{1},alg3par,'UniformOutput',0);

p = blockpanel([commonpar,alg0par,alg1par,alg2par,alg3par]);

% Positive frequency channels
M2 = floor(M/2) + 1;

% Setup blocktream
try
    fs=block(source,varargin{:},'loadind',p,'double');
catch
    % Close the windows if initialization fails
    blockdone(p);
    err = lasterror;
    error(err.message);
end

% Buffer length (30 ms)
bufLen = floor(30e-3*fs);

% Window length in ms
[F,Fdual] = framepair('dgtreal',g,'dual',a,M,'timeinv');
[Fa,Fs] = blockframepairaccel(F,Fdual, bufLen,'segola');

% Prepare buffer
cslicein = zeros(M2,3);
logs = zeros(M2,3);

%%%%% RTISI-LA %%%%%%%
lookback = ceil(M/a) - 1;
framesRTISI = zeros(M,lookback+1+lookaheadmax);
abssRTISI = zeros(M2,1+lookaheadmax);
[gnum,gdnum,specg1,specg2] = comp_rtisilawins(g,gd,a,M);

%%%%% RTPGHI %%%%%%
gamma = -pi/4*gl^2/log(h);
mrange = (1:M2-2)';
tgradmulconst = a*M/gamma/2;
tgradplusconsts = 2*pi*a*mrange/M;
fgradmul = @(fgrad) -gamma/(a*M)*fgrad;
tgradmul = @(tgrad) tgradmulconst*tgrad + tgradplusconsts;

newphase = zeros(M2,2);
tgrad = zeros(M2,3);

%%%%% GSRTISI-LA %%%%%%%%
framesGSRTISI = zeros(M,lookback+1+lookaheadmax);
cframesGSRTISI = zeros(M2,lookback+1+lookaheadmax);
abssGSRTISI = zeros(M2,1+lookaheadmax);
gnums = comp_gsrtisilawins(g,gd,a,M,lookaheadmax);

%%%%% RTPGHI + GSRTISI-LA %%%%%


fola = [];
flag = 1;
%Loop until end of the stream (flag) and until panel is opened
while flag && p.flag

    % Obtain parameters from the control panel
    gain = 10^(p.getParam('GdB')/20); % dB -> val
    do_bypass = p.getParam('Bypass');
    do_causal = 0;
    algno = p.getParam('algno');

    if algno==0
        panelsetvisible(p,alg0parNames,1);
        panelsetvisible(p,alg1parNames,0);
        panelsetvisible(p,alg2parNames,0);
        panelsetvisible(p,alg3parNames,0);
        tol = 10^(-p.getParam('Thr'));
        do_causal = ~p.getParam('rtpghila0');
    elseif algno==1
        panelsetvisible(p,alg0parNames,0);
        panelsetvisible(p,alg1parNames,1);
        panelsetvisible(p,alg2parNames,0);
        panelsetvisible(p,alg3parNames,0);
        lookahead = p.getParam('lookahead1');
        maxit = p.getParam('maxit1');
    elseif algno==2
        panelsetvisible(p,alg0parNames,0);
        panelsetvisible(p,alg1parNames,0);
        panelsetvisible(p,alg2parNames,1);
        panelsetvisible(p,alg3parNames,0);
        lookahead = p.getParam('lookahead2');
        maxit = p.getParam('maxit2');
    elseif algno==3
        panelsetvisible(p,alg0parNames,0);
        panelsetvisible(p,alg1parNames,0);
        panelsetvisible(p,alg2parNames,0);
        panelsetvisible(p,alg3parNames,1);
        lookahead = p.getParam('lookahead3');
        maxit = p.getParam('maxit3');
        do_causal = ~p.getParam('rtpghila3');
    end
    
    % Read block of length bufLen
    [f,flag] = blockread(bufLen);
    % Apply analysis frame
    c = blockana(Fa, f(:,1)*gain);
    
    % This formats the coefficients to M/2+1 x n x W matrix
    cnative = Fa.coef2native(c,size(c));
    
    % Prepare output array of the same size
    cnativeout = zeros(size(cnative));
    
    % Throw away the phase
    s = abs(cnative);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Alter s in any way
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    if do_bypass
        % Bypass
        cnativeout = cnative;
    else
        if algno==0
            % Do RTPGHI
            for ii=1:size(s,2)
                sii= s(:,ii);
                
                % Store abss for later use
                cslicein(:,1:end-1) = cslicein(:,2:end);
                cslicein(:,end) = sii;
                
                siilog = log(sii+eps);
                % siilog(siilog<max(siilog(:))+tt)=tt;
                % Store abss for later use
                logs(:,1:end-1) = logs(:,2:end);
                logs(:,end) = siilog;
                
                % Compute the phase gradient
                idx = 1:2;
                if do_causal
                    fgrad = fgradmul((3*logs(:,3)-4*logs(:,2)+logs(:,1))/2);
                    idx = 2:3;
                else
                    fgrad = fgradmul((logs(:,3)-logs(:,1))/2);
                end
                tgrad(:,1:end-1) = tgrad(:,2:end);
                tgrad(2:end-1,end) = tgradmul(conv2(logs(:,end),[1;0;-1],'valid'));
                
                % Integrate the gradient
                newphase = comp_rtpghiupdate(logs(:,idx),tgrad(:,idx),fgrad,newphase,tol(1),M);  
                cnativeout(:,ii) = sii.*exp(1i*newphase);
            end
        elseif algno==1
            % Do RTISI-LA
            for ii=1:size(s,2)
                sii = s(:,ii);
                framesRTISI(:,1:end-1) = framesRTISI(:,2:end);
                framesRTISI(:,lookback + 1 + lookahead) = 0;
                
                abssRTISI(:,1:end-1) = abssRTISI(:,2:end);
                abssRTISI(:,1 + lookahead) = sii;
                
                [framesRTISI(:,1:lookback + 1 + lookahead),cnativeout(:,ii)] = ...
                    comp_rtisilaupdate(framesRTISI(:,1:lookback + 1 + lookahead),...
                    gnum,specg1,specg2,gdnum,a,M,abssRTISI(:,1:lookahead+1),lookahead,floor(maxit/(lookahead+1)));
            end
        elseif algno==2
            % Do GSRTISI-LA
            for ii=1:size(s,2)
                sii = s(:,ii);
                framesGSRTISI(:,1:end-1) = framesGSRTISI(:,2:end);
                framesGSRTISI(:,lookback + 1 + lookahead) = 0;
                
                cframesGSRTISI(:,1:end-1) = cframesGSRTISI(:,2:end);
                cframesGSRTISI(:,lookback + 1 + lookahead) = 0;
                
                abssGSRTISI(:,1:end-1) = abssGSRTISI(:,2:end);
                abssGSRTISI(:,1 + lookahead) = sii;
                
                [framesGSRTISI(:,1:lookback + 1 + lookahead),...
                 cframesGSRTISI(:,1:lookback + 1 + lookahead),...
                 cnativeout(:,ii)] = ...
                    comp_gsrtisilaupdate(framesGSRTISI(:,1:lookback + 1 + lookahead),...
                                         cframesGSRTISI(:,1:lookback + 1 + lookahead),...
                    gnums,gdnum,a,M,abssGSRTISI(:,1:lookahead+1),lookahead,floor(maxit/(lookahead+1)),0);
            end
        elseif algno==3
            % Do combination of GSRTISI-LA and RTPGHI
                        
            for ii=1:size(s,2)
                sii = s(:,ii);
           
                framesGSRTISI(:,1:end-1) = framesGSRTISI(:,2:end);
                framesGSRTISI(:,lookback + 1 + lookahead) = 0;
                
                cframesGSRTISI(:,1:end-1) = cframesGSRTISI(:,2:end);
                cframesGSRTISI(:,lookback + 1 + lookahead) = 0;
                
                abssGSRTISI(:,1:end-1) = abssGSRTISI(:,2:end);
                abssGSRTISI(:,1 + lookahead) = sii;
                
                cslicein(:,1:end-1) = cslicein(:,2:end);
                cslicein(:,end) = sii;
                
                logs(:,1:end-1) = logs(:,2:end);
                logs(:,end) = log(sii+eps);
               
                % Compute the phase gradient
                idx = 1:2;
                if do_causal || lookahead == 0
                    fgrad = fgradmul((3*logs(:,3)-4*logs(:,2)+logs(:,1))/2);
                    idx = 2:3;
                else
                    fgrad = fgradmul((logs(:,3)-logs(:,1))/2);
                end
                tgrad(:,1:end-1) = tgrad(:,2:end);
                tgrad(2:end-1,end) = tgradmul(conv2(logs(:,end),[1;0;-1],'valid'));
           
               
                if do_causal || lookahead == 0
                    newphase = comp_rtpghiupdate(logs(:,idx),tgrad(:,idx),fgrad,angle(cframesGSRTISI(:,lookback + lookahead)),tol(1),M);  
                    ctmp = abssGSRTISI(:,lookahead+1).*exp(1i*newphase);
                    %ftmp = gdnum.*fftshift(comp_ifftreal(ctmp,M))*M;
                    
                    %framesGSRTISI(:,lookback + 1 + lookahead) = ftmp;
                    cframesGSRTISI(:,lookback + 1 + lookahead) = ctmp;
                    
                    [framesGSRTISI(:,1:lookback + 1 + lookahead),...
                     cframesGSRTISI(:,1:lookback + 1 + lookahead),...
                     cnativeout(:,ii)] = ...
                        comp_gsrtisilaupdate(framesGSRTISI(:,1:lookback + 1 + lookahead),...
                                             cframesGSRTISI(:,1:lookback + 1 + lookahead),...
                        gnums,gdnum,a,M,abssGSRTISI(:,1:lookahead+1),lookahead,floor(maxit/(lookahead+1)),0);
                    
                else
                    newphase = comp_rtpghiupdate(logs(:,idx),tgrad(:,idx),fgrad,angle(cframesGSRTISI(:,lookback + lookahead - 1)),tol(1),M);  
                    ctmp = abssGSRTISI(:,lookahead).*exp(1i*newphase);
                    %ftmp = gdnum.*fftshift(comp_ifftreal(ctmp,M))*M;
                    
                    %framesGSRTISI(:,lookback + lookahead) = ftmp;
                    cframesGSRTISI(:,lookback + lookahead) = ctmp;
                
                    [framesGSRTISI(:,1:lookback + lookahead),...
                     cframesGSRTISI(:,1:lookback + lookahead),...
                     cnativeout(:,ii)] = ...
                        comp_gsrtisilaupdate(framesGSRTISI(:,1:lookback + lookahead),...
                                             cframesGSRTISI(:,1:lookback + lookahead),...
                        gnums,gdnum,a,M,abssGSRTISI(:,1:lookahead),lookahead-1,floor(maxit/(lookahead)),0);
                end
            end            
        end
    end
    
    chat = Fa.native2coef(cnativeout);
    % Apply synthesis frame
    fhat = real(blocksyn(Fs, chat, size(f,1)));
    
    % Analyse again
    %[c2,fola]= blockana(Fa, fhat, fola);
    
    % Plot error
    %blockplot(fobj,Fa,chat);
    % Play the block
    blockplay(fhat);
end
blockdone(p);


function [a,v] = parsevararginfor(what,v,defaultval)
% Parse out a from varargin
apos = find(strcmp(what,v),1);
if ~isempty(apos)
    if apos==numel(v)
        error('%s: Forgotten value for key ''%s''',what,upper(mfilename));
    end
    % Just let ltfatarghelper sort out the correct key-val format
    definput.keyvals.(what) = [];
    [~,~,a] = ltfatarghelper({what},definput,v(apos:apos+1));
    complainif_notposint(a,what,mfilename);
    % And remove it from varargin such that it does not break the rest
    v(apos:apos+1) = [];
else
    a = defaultval;
end


function panelsetvisible(p,params,onoff)

if ischar(params)
    params = {params};
end

cellfun(@(parEl) p.setVisibleParam(parEl,onoff),params,'UniformOutput',0);


