function [f,gnum] = idgtreal_byhand(coef,g,a,M)
% This works like idgtreal (timeinv) but only for FIR windows.

[~,N,W] = size(coef);
L=N*a;
gl = numel(g);
gnum = g;

f = zeros(L,W);
cframes = ifftreal(coef,M)*M;

% This is substituting ifftshift
idxrange = [0:ceil(gl/2)-1,-floor(gl/2):-1];
idxrange2 = mod(idxrange,M) + 1;

for w=1:W
    for n=0:N-1
        idx = mod(n*a + idxrange,L) + 1;
        f(idx,w) = f(idx,w) + cframes(idxrange2,n+1,w).*gnum;
    end
end
