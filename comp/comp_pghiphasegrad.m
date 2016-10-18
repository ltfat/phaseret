function [ tgrad, fgrad, logs ] = comp_pghiphasegrad( s, gamma, a, M, do_timeinv, do_causal )

[M2,N] = size(s);

if do_timeinv
    fgradmul = @(fgrad) -gamma/(a*M)*fgrad;
    tgradmul = @(tgrad) bsxfun(@plus,a*M/gamma*tgrad, 2*pi*a*(0:M2-1)'/M);
else
    fgradmul = @(fgrad) bsxfun(@plus,-gamma/(a*M)*fgrad, -2*pi*a*(0:N-1)/M);
    tgradmul = @(tgrad) a*M/gamma*tgrad;    
end

logs=log(s + eps);
% tt=-10;
% logs(logs<max(logs(:))+tt)=tt;

difforder = 2;


if do_causal
    fgrad = (circshift(logs,[0,2]) -4*circshift(logs,[0,1])+3*logs)/2;
else
    fgrad = pderiv(logs,2,difforder)/size(logs,2);
end

% Undo the scaling done by pderiv and scale properly
tgrad = pderiv(logs,1,difforder)/size(logs,1);

% Fix the first and last rows .. the
% borders are symmetric so the centered difference is 0
tgrad(1,:) = 0;
tgrad(end,:) = 0;

tgrad = tgradmul(tgrad);
fgrad = fgradmul(fgrad);



