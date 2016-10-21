function [partrec,nthframe] = comp_overlayframes(cframes,a,M,n)

N = size(cframes,2);
bufLen = N*a - (a-1) + M-1;
partrec = zeros(bufLen,1);

startidx = floor(M/2);
idxrange = startidx + [0:ceil(M/2)-1,-floor(M/2):-1];
for ii=0:N-1
    idx = ii*a + idxrange + 1;
    partrec(idx) = partrec(idx) + cframes(:,ii+1);
end

if nargin > 3
    nthframe = partrec(a*n+1:a*n + M);
end