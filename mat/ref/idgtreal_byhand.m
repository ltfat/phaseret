function [f,gnum] = idgtreal_byhand(coef,g,a,M)
% This works like idgtreal (timeinv) but only for FIR windows.

[~,N,W] = size(coef);
L=N*a;
gnum = gabwin(g,a,M);

f = zeros(L,W);
cframes = bsxfun(@times,ifftreal(coef,M)*M,gnum);

% This is substituting ifftshift
idxrange = [0:floor(M/2),-ceil(M/2)+1:-1];
for w=1:W
    for n=0:N-1
        idx = mod(n*a + idxrange,L) + 1;
        f(idx,w) = f(idx,w) + cframes(:,n+1,w);
    end
end
