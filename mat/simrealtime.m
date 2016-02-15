function [c,newphase,tgrad,fgrad]=simrealtime(s,lambdaL,a,M,varargin)

thismfilename = upper(mfilename);
complainif_notposint(a,'a',thismfilename);
complainif_notposint(M,'M',thismfilename);

definput.keyvals.tol=[1e-1,1e-10];
definput.keyvals.mask=[];
definput.keyvals.lambda = [];
definput.keyvals.g = [];
definput.flags.phase={'timeinv','freqinv'};
definput.flags.variant={'normal','causal'};
[flags,kv,tol]=ltfatarghelper({'tol','mask'},definput,varargin);
g = kv.g;
[M2,N,W] = size(s);

if W>1
    error('%s: *s* must not be 3 dimensional.',thismfilename);
end

M2true = floor(M/2) + 1;

if M2true ~= M2
    error('%s: Mismatch between *M* and the size of *s*.',thismfilename);
end

L=N*a;

% [~,info]=gabwin(g,a,M,L,'callfun',upper(mfilename));
% 
% if ~info.gauss
%     error(['%s: The window must be a Gaussian window (specified ',...
%            'as a string or as a cell array)'],upper(mfilename));
% end

abss = abs(s);
%lambdaL = info.tfr*L;

if flags.do_timeinv
    fgradmul = @(fgrad) -lambdaL/(a*M)*fgrad;
    tgradmul = @(tgrad) bsxfun(@plus,a*M/lambdaL*tgrad, 2*pi*a*(0:M2-1)'/M);
elseif flags.do_freqinv
    error('%s: FIXME: Not implementet yet.',upper(mfilename));
end

logs=log(abss+realmin);
tt=-10;
logs(logs<max(logs(:))+tt)=tt;

difforder = 2;

if flags.do_normal
    fgrad = pderiv(logs,2,difforder)/size(logs,2);
elseif flags.do_causal
    fgrad = (circshift(logs,[0,2]) -4*circshift(logs,[0,1])+3*logs)/2;
    %fgrad = -1/3*circshift(logs,[0,3]) + 3/2*circshift(logs,[0,2]) -3*circshift(logs,[0,1])+11/6*logs;
    %fgrad = -circshift(logs,[0,1])+logs;
end

% Undo the scaling done by pderiv and scale properly
tgrad = pderiv(logs,1,difforder)/size(logs,1);

% Fix the first and last rows .. the
% borders are symmetric so the centered difference is 0
tgrad(1,:) = 0;
tgrad(end,:) = 0;

tgrad = tgradmul(tgrad);
fgrad = fgradmul(fgrad);
usephase = zeros(M2,2);
tmpmask = zeros(M2,2);
tmpmask(:,1) = 1;
newphase = zeros(M2,N);
c = zeros(M2,N);


for n=1:N
    idx = mod( n-1-1:n-1, N ) + 1;   

    newphasetmp = comp_constructphasereal(abss(:,idx),tgrad(:,idx),fgrad(:,idx),a,M,tol,2,tmpmask,usephase);
    usephase(:,1) = newphasetmp(:,2);
    newphase(:,n) = usephase(:,1);
    
    % Build the coefficients
    c(:,n)=abss(:,n).*exp(1i*newphase(:,n));
end



