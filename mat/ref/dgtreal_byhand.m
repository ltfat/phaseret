function [c,Ls,gnum]=dgtreal_byhand(f,g,a,M)
% This works like dgtreal (timeinv) but only for FIR windows.

[Ls,W] = size(f);
L=dgtlength(Ls,a,a);
N = L/a;
f = postpad(f,L);

gl = numel(g);
gnum = conj(g);

cframes = zeros(M,N,W);

% This is substituting fftshift
idxrange = [0:ceil(gl/2)-1,-floor(gl/2):-1];
idxrange2 = mod(idxrange,M) + 1;

for w=1:W
    for n=0:N-1
        idx = mod(n*a + idxrange,L) + 1;
        cframes(idxrange2,n+1,w) = f(idx,w).*gnum;
    end
end

c = fftreal(cframes);






