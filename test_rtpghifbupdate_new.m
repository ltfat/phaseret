clear all
clc

tol=1e-10;
[f,fs] = gspi;

Ls = length(f);
scales = 2.^(linspace(6,-3.3,200));
[g,a,fc,L, info]=waveletfilters(Ls,scales, 'repeat', 'uniform');

% ufilterbank provides matrix output
corig = ufilterbank(f,g,a);

%If the filterbank is crappy, the result will be crappy. Thus, check
%framebounds!
BoverA = filterbankrealbounds(g,a,L)

% Compute dual filter bank. This is required for reconstruction. 
gd = filterbankrealdual(g,a,L);

% Filterbankconstphase -> Use 'wavelet' flag
[c_constphase,newphase1,usedmask,tgrad,fgrad]=filterbankconstphase(abs(corig), ... %cellfun(@abs,corig,'UniformOutput',0)
    a,info.fc,info.tfr,'tol',1e-10,'wavelet');

% Test filterbankconstphase result
Emags = magnitudeerrdb(corig,c_constphase) % Foolproofing

f_rec_constphase = ifilterbank(c_constphase,gd,a,'real');
c_rec_constphase = ufilterbank(f_rec_constphase,g,a);
Emags_constphase = magnitudeerrdb(corig,c_rec_constphase) % Spectral error (below -20 is good!)

%generate same spectrum as used in filterbankconstphase (with "natural"
%scaling)
[N,M] = size(corig);

asan = comp_filterbank_a(a,M);
a = asan(:,1)./asan(:,2);
scal = 1./(sqrt(asan(:,1)./asan(:,2)));
abss = abs(corig.').*repmat(scal, 1, N);

%transposition as of yet still necessary (will change that later)
tgrad2 = tgrad.';
fgrad2 = fgrad.';

newphase = zeros(M,N);
crec = zeros(M,N);

%scale fgrad and tgrad to match the normalized fc
fgrad2 = -pi*fgrad2;
tgrad2 = a(1)*(tgrad2 + info.fc') * pi;%scale to samples

for n=1:N
    idx = mod( n-1-1:n-1, N ) + 1;
    nprev = mod( n-2, N ) + 1;
    %logarithm could be a problem
    newphase(:,n) = comp_rtpghifbupdate(abss(:,idx),info.fc,tgrad2(:, idx), fgrad2(:,n), newphase(:,nprev),tol,numel(fc));%(numel(fc)-1)*2);
    crec(:,n)=abss(:,n).*exp(1i*newphase(:,n));
end

% Test results
E = magnitudeerrdb(corig,crec.'); % Wrong!

% Undo scaling? Yes!!!
%crec = (crec./repmat(scal, 1, N)).';
crec = crec.';
%Emags = magnitudeerrdb(corig,crec) % Foolproofing

f_rec_rtpghi = ifilterbank(crec,gd,a,'real');
c_rec_rtpghi = ufilterbank(f_rec_constphase,g,a);
Emags_rtpghi = magnitudeerrdb(corig,c_rec_rtpghi) % Spectral error (below -20 is good!)

%% Uncomment to listen
% soundsc(f_rec_rtpghi, fs)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
% soundsc(f_rec_constphase, fs)

%a few further debugging things...

% Did this above already!
% gd=filterbankdual(g,a,L); %dual makes it worse?
% newphase2 = reshape(cell2mat(newphase1), length(newphase1{1}), length(cell2mat(newphase1))/length((newphase1{1})));
% 
% 
% figure(1)
% subplot(2,1,1)
% plot(sum(newphase))
% xlabel('rtpghifb')
% hold on
% subplot(2,1,2)
% plot(sum(newphase2,2))
% xlabel('filterbankconstphase')
% 
%cprojconst = filterbank(ifilterbankiter(c,g,a,'pcg','tol',1e-6),g,a);
%Cdbfbconst = 20*log10( norm(abs(cell2mat(corig)) - abs(cell2mat(cprojconst)) )/norm( abs(cell2mat(corig))) );
% 
%cproj = filterbank(ifilterbankiter(crec.',g,a,'pcg','tol',1e-6),g,a);
%Cdb = 20*log10( norm(abs(cell2mat(corig)) - abs(cell2mat(cproj)) )/norm( abs(cell2mat(corig))) );


   
