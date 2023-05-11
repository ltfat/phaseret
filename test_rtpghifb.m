clear all
clc

tol=1e-10;
[f,fs] = gspi;

Ls = length(f);
scales = 2.^(linspace(6,-3.3,200));
[g,a,fc,L, info]=waveletfilters(Ls,scales, 'repeat', 'uniform');
corig = ufilterbank(f,g,a);

[c,newphase,tgrad,fgrad]=rtpghifbwl(abs(corig.'),a(1),info.fc,info.tfr);

c = c.';
Emags = magnitudeerrdb(corig,c) % Foolproofing

gd = filterbankrealdual(g,a,L);
f_rec_rtpghi = ifilterbank(c,gd,a,'real');
c_rec_rtpghi = ufilterbank(f_rec_rtpghi,g,a);
Emags_rtpghi = magnitudeerrdb(corig,c_rec_rtpghi)

