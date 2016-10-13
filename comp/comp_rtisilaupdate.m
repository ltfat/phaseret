function [cframes2,c] = comp_rtisilaupdate(cframes,gnum,specg1,specg2,gdnum,a,M,sframes,lookahead,maxit)
% cframes: lookback + 1 + lookahead == size(cframes,2)
% cframes2: 1 + lookahead == size(cframes2,2)
% sframes: 1 + lookahead == size(sframes,2)

N = size(cframes,2);
lookback = N - lookahead - 1;


%% 2) Other lookahead frames and the submit frame
% for nback = lookahead:-1:0
%     indx = lookback+nback+1;
%     
%     if nback == lookahead
%         prd = overlaynthframe(cframes,indx,specg1,a,M);
%     else
%         prd = overlaynthframe(cframes,indx,gnum,a,M);
%     end
%     
%     [cframes(:,indx),c] = phaseupdate(prd,sframes(:,nback+1),gdnum,M);
% end

frameorder = lookahead:-1:0;
%% 3) All the other iterations
for iter=1:maxit
    for nback = frameorder
        indx = lookback+nback+1;
        
        if nback == lookahead
            if iter == 1
                prd = overlaynthframe(cframes,indx,specg1,a,M);
            else
                prd = overlaynthframe(cframes,indx,specg2,a,M);
            end
        else
            prd = overlaynthframe(cframes,indx,gnum,a,M);
        end
        [cframes(:,indx),c] = phaseupdate(prd,sframes(:,nback+1),gdnum,M);
    end
end

% modd = ones(size(c));
% modd(2:2:end) = -1;
% c = c.*modd;
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

c = s.*exp(1i*angle(comp_fftreal(circshift(cbuf,-floor(M/2)))));
frame = gdnum.*(circshift(comp_ifftreal(c,M),floor(M/2)))*M;
