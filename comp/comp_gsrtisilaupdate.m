function [cframes2,coefbuf,c] = comp_gsrtisilaupdate(cframes,coefbuf,gnums,gdnum,a,M,sframes,lookahead,maxit,do_energyorder)
% cframes: lookback + 1 + lookahead == size(cframes,2)
% cframes2: 1 + lookahead == size(cframes2,2)
% sframes: 1 + lookahead == size(sframes,2)

N = size(cframes,2);
M2 = floor(M/2) + 1;
lookback = N - lookahead - 1;

if do_energyorder
    % New frame order, force the commit frame update to be the last
    [~,frameorder] = sort(sum(sframes(:,2:end)),'descend');
    frameorder = [frameorder, 0];
else
    frameorder = lookahead:-1:0;
end

%% 3) All the other iterations
for iter=1:maxit
    for nback = frameorder
        indx = lookback+nback+1;

        prd = overlaynthframe(cframes,indx,gnums(:,nback+1),a,M);
        [cframes(:,indx), coefbuf(:,indx)] = phaseupdate(prd,sframes(:,nback+1),gdnum,M);
        %c = coefbuf(:,indx);
    end
end

% ifftshift 
modd = ones(M2,lookahead + 1);
modd(2:2:end,:) = -1;
coefbuf(:,lookback + 1:end) = coefbuf(:,lookback + 1:end).*modd;
c = coefbuf(:,lookback + 1);

cframes2 = cframes;


function partrec = overlaynthframe(cframes,n,gnum,a,M)
% cframes are reconstructed frames right before overlap add

ceilM2 = ceil(M/2);
floorM2 = floor(M/2);
%cframes = cframes([floorM2+1:M, 1:floorM2],:);

N = size(cframes,2);
% Initial sum, just the n-th frame itself
partrec = cframes(:,n);

% Go trough frames to the right
for ii = n+1:N
    jj = (ii - n)*a + 1;
    nSamp = M - jj + 1;
    if nSamp <= 0
        break;
    end
    partrec(jj:end) = partrec(jj:end) + cframes(1:nSamp,ii);
end

% Go trough frames to the left
for ii = n-1:-1:1
    jj = (n-ii)*a + 1;
    nSamp = M - jj+1;
    if nSamp <= 0
        break;
    end
    partrec(1:nSamp) = partrec(1:nSamp) + cframes(jj:end,ii);
end

%partrec = fftshift(partrec,1);
%partrec = partrec([ceilM2+1:M, 1:ceilM2]);

if ~isempty(gnum)
    partrec = partrec.*gnum;
end

function [frame,c] = phaseupdate(cbuf,s,gdnum,M)

c = s.*exp(1i*angle(comp_fftreal((cbuf))));
frame = gdnum.*(comp_ifftreal(c,M))*M;
