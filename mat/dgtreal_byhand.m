function [c,Ls,gnum]=dgtreal_byhand(f,g,a,M)
[Ls,W] = size(f);
L=dgtlength(Ls,a,M);
N = L/a;
f = postpad(f,L);

gnum = conj(gabwin(g,a,M));
cframes = zeros(M,N,W);

% This is substituting fftshift
idxrange = [0:floor(M/2),-ceil(M/2)+1:-1];

for w=1:W
    for n=0:N-1
        idx = mod(n*a + idxrange,L) + 1;
        cframes(:,n+1,w) = f(idx,w);
    end
end

c = fftreal(bsxfun(@times,cframes,gnum));






