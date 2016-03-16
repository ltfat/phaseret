function demo_blockproc_phaseret(source,varargin) %RUNASSCRIPT
%DEMO_BLOCKPROC_PHASERET STFT phase reconstruction algorithms comarison
%   Usage: demo_blockproc_phaseret('gspi.wav')
%          demo_blockproc_phaseret(...,'a',a)` 
%
%   For additional help call |demo_blockproc_phaseret| without arguments.
%
%   This file does real-time reconstruction from the STFT magnitude using
%   3 algorithms. No effect, just plain reconstruction.
%       
%   `demo_blockproc_phaseret(...,'a',a)` allows setting the hop factor.
%   Default value is 256.
%
%   .. image:: ../images/phaseret_rtpghi.png
%   .. image:: ../images/phaseret_spsi.png
%   .. image:: ../images/phaseret_rtisi.png
% 
%   References: ltfatnote043
%   References: ltfatnote043

% AUTHOR: Zdenek Prusa

if demo_blockproc_header(mfilename,nargin)
    return;
end

%fobj = blockfigure();

% Parse out a from varargin
apos = find(strcmp('a',varargin),1);
if ~isempty(apos)
    if apos==numel(varargin)
        error('%s: Forgotten value for key ''a''',upper(mfilename));
    end
    % Just let ltfatarghelper sort out the correct key-val format
    definput.keyvals.a = [];
    [~,~,a] = ltfatarghelper({'a'},definput,varargin(apos:apos+1));
    complainif_notposint(a,'a',mfilename);
    % And remove it from varargin such that it does not break the rest
    varargin(apos:apos+1) = [];
else
    a = 256;
end

% Number of frequency channels
M = 2048;
gl = 2048;
% Window support is M
g = {'blackman',gl};

lookaheadmax = 16;

% Control pannel (Java object)
% Each entry determines one parameter to be changed during the main loop
% execution.
commonpar = {
    {'GdB','Gain',-20,20,0,21},...
    {'Bypass' ,'Bypass' ,0,1,0,2},...
    {'algno' ,'Algorithm' ,0,2,0,3},...
    };
... % Alg 0 params:
    alg0par = {
    {'Thr','Thr',0,16,6,100},...
    {'causal' ,'Causal' ,0,1,0,2},...
    };
alg0parNames = cellfun(@(aEl) aEl{1},alg0par,'UniformOutput',0);

alg2par = {
    {'lookahead','Lookahead',0,lookaheadmax,0,lookaheadmax+1},...
    {'maxit','Max. iter.',1,16,8,16},...
    };
alg2parNames = cellfun(@(aEl) aEl{1},alg2par,'UniformOutput',0);

p = blockpanel([commonpar,alg0par,alg2par]);



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

lookback = ceil(M/a) - 1;
framesRTISI = zeros(M,lookback+1+lookaheadmax);
abssRTISI = zeros(M2,1+lookaheadmax);

gnum = normalize(gabwin(g,a,M),'inf');
% Synthesis window
gd = gabdual(g,a,M);
% ... as array
gdnum = gabwin(gd,a,M);
wins = repmat(gdnum,1,2*lookback+1);
specg1 = overlayframes(wins,a,M);
specg1 = (specg1(1:M));
wins(:,1) = 0;
specg2 = overlayframes(wins,a,M);
specg2 = (specg2(1:M));
gnum = fftshift(gnum);
gdnum = fftshift(gdnum);

lambdaL = gl^2*0.17937;
%lambdaL = -pi*gl^2/(4*log(h));
mrange = (1:M2-2)';
fgradmul = @(fgrad) -lambdaL/(a*M)*fgrad;
tgradmul = @(tgrad) bsxfun(@plus,a*M/lambdaL*tgrad, 2*pi*a*mrange/M);

tmpmask = zeros(M2,2);
tmpmask(:,1) = 1;
newphase = zeros(M2,2);

fgrad = zeros(M2,2);
tgrad = zeros(M2,2);

startphase = zeros(M2,1);
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
        panelsetvisible(p,alg2parNames,0);
        tol = 10^(-p.getParam('Thr'));
        do_causal = p.getParam('causal');
    elseif algno==1
        panelsetvisible(p,alg0parNames,0);
        panelsetvisible(p,alg2parNames,0);
    elseif algno==2
        panelsetvisible(p,alg0parNames,0);
        panelsetvisible(p,alg2parNames,1);
        lookahead = p.getParam('lookahead');
        maxit = p.getParam('maxit');
    end
    
    % Read block of length bufLen
    [f,flag] = blockread(bufLen);
    % Apply analysis frame
    c = blockana(Fa, f*gain);
    
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
            % Do heap integration
            for ii=1:size(s,2)
                sii= s(:,ii);
                
                % Store abss for later use
                cslicein(:,1:end-1) = cslicein(:,2:end);
                cslicein(:,end) = sii;
                
                siilog = log(sii+realmin);
                % siilog(siilog<max(siilog(:))+tt)=tt;
                % Store abss for later use
                logs(:,1:end-1) = logs(:,2:end);
                logs(:,end) = siilog;
                
                % Compute the phase gradient
                idx = 1:2;
                if do_causal
                    fgrad(:,2) = fgradmul((3*logs(:,3)-4*logs(:,2)+logs(:,1))/2);
                    idx = 2:3;
                else
                    fgrad(:,2) = fgradmul((logs(:,3)-logs(:,1))/2);
                end
                tgrad(2:end-1,:) = tgradmul(conv2(logs(:,idx),[1;0;-1],'valid')/2);
                
                % Integrate the gradient
                tmp = comp_constructphasereal(cslicein(:,idx),tgrad,fgrad,a,M,tol(1),2,tmpmask,newphase);
                newphase(:,1) = tmp(:,2);
                
                % Use the new phase
                cnativeout(:,ii) = sii.*exp(1i*tmp(:,2));
            end
        elseif algno==1
            % Do SPSI
            [cnativeout,startphase] = comp_spsi(s,a,M,startphase);
        elseif algno==2
            % Do RTISI-LA
            for ii=1:size(s,2)
                sii = s(:,ii);
                framesRTISI(:,1:end-1) = framesRTISI(:,2:end);
                framesRTISI(:,lookback + 1 + lookahead) = 0;
                
                abssRTISI(:,1:end-1) = abssRTISI(:,2:end);
                abssRTISI(:,1 + lookahead) = sii;
                
                [framesRTISI(:,1:lookback + 1 + lookahead),cnativeout(:,ii)] = ...
                    comp_rtisilaupdate(framesRTISI(:,1:lookback + 1 + lookahead),...
                    gnum,specg1,specg2,gdnum,a,M,abssRTISI(:,1:lookahead+1),lookahead,maxit);
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



function panelsetvisible(p,params,onoff)

if ischar(params)
    params = {params};
end

cellfun(@(parEl) p.setVisibleParam(parEl,onoff),params);

function partrec = overlayframes(cframes,a,M)

N = size(cframes,2);
bufLen = N*a - (a-1) + M-1;
partrec = zeros(bufLen,1);

startidx = ceil(M/2)-1;
idxrange = startidx + [0:floor(M/2),-ceil(M/2)+1:-1];
for n=0:N-1
    idx = n*a + idxrange + 1;
    partrec(idx) = partrec(idx) + cframes(:,n+1);
end
