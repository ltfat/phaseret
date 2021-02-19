%DEMO_PHASEVOCODER Vocoder-based time-scale and pitch modification 
%
%   In this file, we demonstrate how to use |phasevocoder| to perform 
%   time-scale and pitch modification of audio. 
%   
%   This script is comprised of 5 parts and which can separately be adapted 
%   by the user: 
%
%   1) An audio signal is loaded and (by default) truncated to 2 seconds length. 
%   2) Set Vocoder-related parameters, in particular stretching factors, 
%      used algorithms and STFT parameters.
%   3) Perform time-scale and pitch modification [only modify this if you
%      are an experienced user].
%   4) Display instructions for accessing the results [no need to change
%      this at all].
%   5) Plot waveform and DGT spectrogram of the input and one example from 
%      the pitch modified signals.
%
%   .. figure:: 
%
%      Input signal and pitch-shifted signal 
%
%      Waveform and DGT spectrogram of the input signal (top) and a
%      pitch-shifted signal (bottom). 
%
%   See also: phasevocoder, rtpghi
%
%   References: ltfatnote050, lado99b

% AUTHORS: Zdenek Prusa, Nicki Holighaus and Clara Hollomey

disp('Type "help demo_phasevocoder" to see a description of how this demo works.');

%% Load and prepare input signal
[insig,fs]=gspi; % Load input signal 
Ls= min(size(insig,1),2*fs); % Restrict signal to at most 2 seconds
insig = sum(insig(1:Ls,:),2)/size(insig,2); % Merge signal channels

%% Vocoder-related parameters
stretches = [1/4,1/3,1/2,2/3,3/2,2,3,4]; % Desired stretch factors
algorithms = {'pv','rtpghi'};

% For vocoders using the full phase gradient, it is recommended that 
% M >= max(gl,max(stretch*gl))  
M = 8192;  % Number of STFT channels
gl = 2048; % STFT window length
wintype = 'hann'; % Window type, see firwin for available types

% Set analysis and synthesis hop size. For stretch factors > 1, a_syn will 
% be used, for factors below 1, a_ana  will be used. 
a_syn = 256;
a_ana = 256;

% Refer to the help of 'phasevocoder' for more information on the
% vocoder-related parameters.

%% Initialization and signal modification

% Initialize cell array to store time-scale modified signals
outsigs_stretched = cell(numel(stretches),numel(algorithms));
% Initialize cell array to store pitch-shifted signals
outsigs_shifted = cell(numel(stretches),numel(algorithms));

% Compute time-scale and pitch modified signals
for kk = 1:numel(stretches)
    fprintf('\n')
    disp(['Computing time-scale and pitch modified signals for stretch factor ',...
         num2str(kk), ' of ',num2str(numel(stretches)),': ',num2str(stretches(kk),'%0.3f')]);
    for jj = 1:numel(algorithms)
        disp(['Doing so for algorihms ',...
         num2str(jj), ' of ',num2str(numel(algorithms)),': ',algorithms{jj}]);
        if stretches(kk) >= 1
            outsigs_stretched{kk,jj} = phasevocoder(insig,stretches(kk),gl,M,wintype,'a_syn',a_syn,algorithms{jj});
            % A pitch-shifted signal can be computed by resampling time-scale
            % modified signal to obtain a signal shifted by
            % 12*log2(stretch) semitones:
            %
            % outsigs_shifted{kk,jj} = dctresample(outsigs_stretched{kk,jj},Ls); 
            %
            % or directly by calling `phasevocoder` with the 'shift' flag:
            
            outsigs_shifted{kk,jj} = phasevocoder(insig,12*log2(stretches(kk)),gl,M,'shift',wintype,'a_syn',a_syn,algorithms{jj});
        else
            outsigs_stretched{kk,jj} = phasevocoder(insig,stretches(kk),gl,M,wintype,'a_ana',a_ana,algorithms{jj});
            %
            % outsigs_shifted{kk,jj} = dctresample(outsigs_stretched{kk,jj},Ls);    
            %
            outsigs_shifted{kk,jj} = phasevocoder(insig,12*log2(stretches(kk)),gl,M,'shift',wintype,'a_ana',a_ana,algorithms{jj});
        end
    end
end

%% Display information for accessing the results

% Preparation
algstring = [];
for jj = 1:numel(algorithms)
    algstring = strcat(algstring,[' ',num2str(jj),': ',algorithms{jj},', ']);
end
algstring = algstring(1:end-1);

stretchstring = [];
pitchstring = [];
for kk = 1:numel(stretches)
    stretchstring = strcat(stretchstring,[' ',num2str(kk),': ',num2str(stretches(kk),'%0.3f'),',']);
    pitchstring = strcat(pitchstring,[' ',num2str(kk),': ',num2str(12*log2(stretches(kk)),'%0.3f'),',']);
end
stretchstring = stretchstring(1:end-1);
pitchstring = pitchstring(1:end-1);

% Actually display instructions
fprintf('\n')
disp('Time-scale modified signals are stored in: outsigs_stretched');
disp(['Columns correspond to algorithm:',algstring]);
disp(['Rows correspond to stretch factor:',stretchstring]);
fprintf('\n')
disp('For playback of a modified signal run: soundsc(outsigs_stretched{k,j},fs)');
disp('For playback of the input run: soundsc(insig,fs)');
fprintf('\n')

fprintf('\n')
disp('Pitch modified signals are stored in: outsigs_shifted');
disp(['Columns correspond to algorithm:',algstring]);
disp(['Rows correspond to pitch-shift in semitones:',pitchstring]);
fprintf('\n')
disp('For playback of a modified signal run: soundsc(outsigs_shifted{k,j},fs)');
disp('For playback of the input run: soundsc(insig,fs)');
fprintf('\n')

%% Plot example

figure(1);
subplot(2,2,1);
plot(insig);
ylim([-1 1]);
xlabel('Time in samples');
ylabel('Signal amplitude');
title('Input signal');

subplot(2,2,2);
c = dgtreal(insig,{wintype,gl},a_ana,M);
plotdgtreal(c,a_ana,M,fs,60);
title('DGT spectrogram of input signal');

subplot(2,2,3);
plot(outsigs_shifted{end,end});
ylim([-1 1]);
xlabel('Time in samples');
ylabel('Signal amplitude');
title(['Pitch modified signal shifted by ',num2str(12*log2(stretches(end)),'%0.3f'),...
' semitones with algorithm: ',algorithms{2}]);

subplot(2,2,4);
c = dgtreal(outsigs_shifted{end,end},{wintype,gl},a_ana,M);
plotdgtreal(c,a_ana,M,fs,60);
title('DGT spectrogram of pitch modified signal');
