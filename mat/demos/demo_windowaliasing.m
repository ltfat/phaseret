a = 512;
M = 2048;
gl = M;
g = {'hann',gl};
L = 10*M;

b = L/M;
N = L/a;

gnum = normalize(gabwin(g,a,M),'inf');
gnum2 = normalize(fir2long(gnum,L),'inf');
gnumFFT2 = normalize(fft(gnum2),'inf'); 

gnumAl = sum(reshape(gnum2,M,b),2);
gnumFFTAl = sum(reshape(gnumFFT2,N,a),2);

figure(1);
subplot(2,1,1);plot(fftshift(real([gnum,gnumAl])));axis tight;
gnumAlShort = long2fir(gnumFFT2,N);

subplot(2,1,2);plot(fftshift(real([gnumAlShort,gnumFFTAl])));axis tight;