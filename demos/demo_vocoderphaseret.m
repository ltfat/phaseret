%DEMO_VOCODERPHASERET Vocoder-based time stretching and picth shifting 
%
%   The scripts shows vocoder ...
%
%   .. figure:: 
%
%      Title
%
%      Description
%
%   Please note that if only pitch-shifting is resired, it is advantenous
%   to do resampling before time stretching if pitch shifting upwards.
%
%   References: lado99b
%

% AUTHOR: Zdenek Prusa
%

% Load a signal
[f,fs]= gspi; Ls= numel(f);
f = f(:,1);

% Choose either pitch shift in semitones (even nonintegers) ...
semitoneshift = -6;
% ... or the time scale ratio directly
timescalerat = 1/(2^(semitoneshift/(12)));

% Number of frequency channels
M = 2048;
% Analysis hop factor
a = 256;
% Window used
g = {'blackman',M};
% Synthesis hop factor
newa = floor23(a/timescalerat);

% Determine compatible length
Lsmallestnewa = dgtlength(1,newa,M);
Lsmallesta = dgtlength(1,a,M);
Lsmallest = dgtlength(1,Lsmallestnewa,Lsmallesta);
L = dgtlength(Ls,Lsmallest,M);
f = postpad(f,L);

% Analysis
[c,~,gnum] = dgtreal(f,g,a,M,'timeinv');

% Adjust to compatible size
N = size(c,2);
newN = ceil(N*newa/Lsmallestnewa)*Lsmallestnewa/newa;
c = postpad(c,newN,0,2);

% Phase reconstruction
gl = numel(gnum);
[chatint] = pghi(abs(c),0.17954*gl^2,newa,M,'timeinv');

% Synthesis
fscale = idgtreal(chatint,{'dual',g},newa,M,'timeinv');
fscale = fscale(1:floor(numel(f)*newa/a));

% Resample to the original duration (more or less)
fshift = dctresample(fscale,newN*newa);


disp('To play the original run: soundsc(f,fs)');
disp('To play the time stretched/compressed version run: soundsc(fscale,fs)');
disp('To play the pitch shifted version run: soundsc(fshift,fs)');



