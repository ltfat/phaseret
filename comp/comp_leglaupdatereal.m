function [cprojupd,cproj] = comp_leglaupdaterealm(c,kern,s,a,M,do_onthefly)
% Mods 0 - just 

if nargin < 6
    do_onthefly = 0;
end

[M2,N] = size(c);

kern = ifftshift(kern);


if rem(M,2)==1
   cfull = [c;conj(c(end:-1:2,:))];
else
   cfull = [c;conj(c(end-1:-1:2,:))];
end

cprojupd = c;
cproj = zeros(size(c));
[kernh,kernw] = size(kern);


if do_onthefly
    if 1
    for n=1:N
        kernn = phasekernfi(kern,(n-1),a,M);
        %kernn = kern;
        idxn = mod(n-1 + fftindex(kernw),N)+1;

        for m=1:M2
            idxm = mod(m-1 + fftindex(kernh),M)+1;
            % Project a single coefficient
            cproj(m,n) = sum(sum(kernn.*cfull(idxm,idxn)));
            % Use only the angle from the projection
            cprojupd(m,n) = s(m,n)*exp(1i*angle(cproj(m,n)));
            % Update the input so that it is used by neighbors
            cfull(m,n) = cprojupd(m,n);
        end
    end
    else
        
        [ssort,idx] = sort(s(:),'descend');
        
        %idx = idx(ssort<1e-2*ssort(1));
        colIds = floor((idx-1)/M2)+1;
        rowIds = mod((idx-1),M2)+1;
        for ii = 1:numel(idx)
            n = colIds(ii); m = rowIds(ii);
            kernn = phasekernfi(kern,(n-1),a,M);
            %kernn = kern;
            idxn = mod(n-1 + fftindex(kernw),N)+1;
            idxm = mod(m-1 + fftindex(kernh),M)+1;
            % Project a single coefficient
            cproj(m,n) = sum(sum(kernn.*cfull(idxm,idxn)));
            % Use only the angle from the projection
            cprojupd(m,n) = s(m,n)*exp(1i*angle(cproj(m,n)));
            % Update the input so that it is used by neighbors
            cfull(m,n) = cprojupd(m,n);
        end
    end
else
    for n=1:N
        kernn = phasekernfi(kern,(n-1),a,M);
        %kernn = kern;
        idxn = mod(n-1 + fftindex(kernw),N)+1;

        for m=1:M2
            idxm = mod(m-1 + fftindex(kernh),M)+1;
            cproj(m,n) = sum(sum(kernn.*cfull(idxm,idxn)));
        end    
    end
    cprojupd = s.*exp(1i*angle(cproj(1:M2,:)));
end



function kernm = phasekernfi(kern,n,a,M)

kernh = size(kern,1);
l = 2*pi*n*fftindex(kernh)*a/M;

kernm = bsxfun(@times,kern,exp(1i*l));

% M/a periodic in m
function kernm = phasekernti(kern,m,a,M)

[kernh, kernw] = size(kern);
l = 2*pi*m*fftindex(kernw)'*a/M;

kernm = bsxfun(@times,kern,exp(1i*l));
