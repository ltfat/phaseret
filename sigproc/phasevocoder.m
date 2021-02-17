function [fout, info] = phasevocoder(f,stretch,varargin)
%PHASEVOCODER  Vocoder-based time-stretching and pitch-shifting 
%   Usage:  fout=phasevocoder(f,stretch);
%           fout=phasevocoder(f,stretch,...);
%           [fout, info]=phasevocoder(...);
%
%   Input parameters:
%         f        : Input signal.
%         stretch  : Desired scalar stretch factor.
%   Output parameters:
%         fout     : Time-stretched/time-compressed output signal.
%         info     : Meta-data - Parameters and STFT analysis/synthesis
%                    coefficients
%
%   `fout=phasevocoder(f,stretch)` computes from the input signal *f* a
%   time-stretched (`stretch > 1`) or time-compressed (`stretch<1`) output 
%   signal *fout* of length `stretch*length(f)`. In contrast to simple 
%   resampling, vocoder-based time-stretching attempts to preserve the 
%   frequency characteristics of the input. The latter is achieved by phase
%   manipulation. The stretch factor must be positive. 
%   
%   If the additional flag `shift` is specified, pitch.shifting is 
%   performed instead, e.g., `fout=phasevocoder(f,stretch,'shift')`. Here, 
%   *stretch* is interpreted as the desired pitch-shift in semitones. 
%   Pitch-shifting is obtained by time-stretching the input with factor 
%   `2^{stretch/12}`, followed by re-sampling to the input signal length. 
%   Resampling is achieved by using |dctresample|. 
%
%   The second output argument *info* is a struct that contains various
%   meta-data related to the internally used short-time Fourier transform: 
%   The analysis and synthesis hop sizes *a_ana*, *a_syn*, the number of 
%   frequency channels *M*, the analysis coefficients *c_ana* of the input 
%   signal *f* and the synthesis coefficients *c_syn* of the time-scale 
%   modified output signal *fout*. If pitch-shifting is performed, *c_syn*
%   corresponds to the output before resampling.
%
%   The desired window length *gl* and number of channels *M* for the 
%   required short-time Fourier transform can optinally be passed as
%   position-dependent parameters:
%
%      fout=phasevocoder(f,stretch,gl);
%      fout=phasevocoder(f,stretch,gl,M);
% 
%   See below for more information.
%
%   Additional parameters
%   ---------------------
%
%   Currently, `phasevocoder` supports two phase vocoder implementations,
%   selected by the following flags:
%
%      'rtpghi'   The 'Phase Vocoder Done Right' proposed by Prusa and 
%                 Holighaus. This vocoder scales the phase gradient
%                 according to the stretch  and uses real-time phase
%                 gradient heap integration to construct the output phase, 
%                 see  |rtpghi| for details. This is the default.
%
%      'pv'       The classic phase vocoder, originally proposed by 
%                 Portnoff. In the classic phase vocoder, only the
%                 time-direction phase derivative is scales according to
%                 the stretch factor. The output phase is obtained by
%                 channel-wise trapezoidal integration of the scaled phase 
%                 derivative.
%
%   By providing the optional argument 'stretch' (default) or 'shift', the
%   user can switch between time-scale modification ('stretch') and pitch
%   manipulation ('shift').
%
%   `phasevocoder` allows to specify the window type used in the internal 
%   short-time Fourier transform by providing any window type supported by 
%   |firwin| as optional flag. See |firwin| for available types. A 'hann'
%   window is used by default.
%
%   `phasevocoder` allows to specify the manner in which the required phase 
%   derivatives are obtained. Options: 'phase', 'dgt', 'abs'. See
%   |gabphasederivs| for more information. When the 'phase' method is
%   selected, it is possible to specify the used finite difference scheme:
%   'forward', 'backward' (default), or 'centered'.
% 
%   More implementations are planned for future releases.
%
%   `phasevocoder` accepts the following optional parameters:
%
%     'gl',gl            Length of the analysis window; the default is 
%                        *gl=4096*. 
%
%     'M',M              Number of frequency channels. In 'rtpghi' mode and 
%                        for moderate stretch factors, it is recommended 
%                        that `M >= max(gl,stretch*gl)` for optimal 
%                        synthesis quality; the default is *M=2048*.
% 
%     'a_ana',a_ana      Analysis hop size; derived from synthesis hop size
%                        *a_syn* and stretching factor by default. 
%
%     'a_syn',a_syn      Synthesis hop size; the default is *a_syn=512*.
%                        This value is ignored if *a_ana* is specified. 
%                        In that case, *a_syn* is derived from *a_ana* and 
%                        *stretch* instead.  
%
%   See also:  rtpghi, firwin, dgtreal
%
%   References: ltfatnote050, ltfatnote048, po78, lado99b
%

% AUTHOR: Zdenek Prusa, Nicki Holighaus, Clara Hollomey

definput.import = {'firwin'};
definput.keyvals.M = 2048;
definput.keyvals.a_ana = [];
definput.keyvals.a_syn = 512;
definput.keyvals.gl = 4096;
definput.flags.method = {'rtpghi','pv'};
definput.flags.gradmethod = {'phase','dgt','abs'};
definput.flags.diffmethod = {'backward','centered', 'forward'};
definput.flags.shift = {'stretch','shift'};
definput.keyvals.dim = [];

[flags,kv,gl,M]=ltfatarghelper({'gl','M'},definput,varargin);
[f,~,Ls]=assert_sigreshape_pre(f,[],kv.dim,upper(mfilename));

if flags.do_shift
    stretch = 2^(stretch/12);
end

if stretch <= 0 
    error(['%s: Stretch factor must be larger than zero. Stretch factor given: %0.3g'],...
              upper(mfilename),stretch);
elseif stretch > 8 || stretch < 1/8
    warning(['This implementation is intended for use with stretch factors '...
        'between 0.125 and 8. Stretch factor given: %0.3g'],stretch); 
end

if ~isempty(kv.a_ana)
    a_ana = kv.a_ana;
    a_syn = floor(a_ana*stretch);
else
    a_syn = kv.a_syn;
    a_ana = floor(a_syn/stretch);
end

% Check that neither analysis nor synthesis hop sizes are larger than
% the window length. Violation leads to an incomplete DGT.
if a_ana >= gl
    error(['%s: Analysis hop size a_ana = %d is larger than window length ',...
          'gl = %d. Please increase window length manually.'],...
              upper(mfilename),a_ana,gl);
elseif a_syn >= gl
    error(['%s: Synthesis hop size a_syn = %d is larger than window length ',...
          'gl = %d. Please increase window length manually.'],...
              upper(mfilename),a_syn,gl);
end

% Determine suitable input signal length. Must be a multiple of a_ana. To
% prevent periodic border artifacts, the window length is accounted for.
% Using DGTLENGTH in this way circumvents the usual requirement that L is
% an integer multiple of a_ana and M.
if flags.do_shift && stretch > 1 % If pitch-shifting upwards, resample first
    f = dctresample(f,floor(Ls/stretch));
    Ls = floor(Ls/stretch);
end
L = dgtlength(max([Ls+gl,kv.M]),a_ana,a_ana);

b_ana = L/kv.M; % Frequency hop size
N = L/a_ana; % Number of time steps
M2 = floor(kv.M/2) + 1; % Actual number of channels for DGTREAL
%truestretch = a_syn/a_ana;

% Determine length of stretched signal
exactnewL = floor(Ls*stretch);

% Determine suitable synthesis length. Must be a multiple of a_syn. To
% prevent periodic border artifacts, the window length is accounted for.
% Using DGTLENGTH in this way circumvents the usual requirement that L is
% an integer multiple of a_ana and M.
newL = dgtlength(max([Ls*stretch+gl,kv.M]),a_syn,a_syn);
newN = newL/a_syn;

% Extend f to length L 
f = postpad(f,L);
% Initialize Array of processed coefficients
chat = [];

% Initialize analysis and synthesis windows
g = {flags.wintype,gl,'inf'};
g_ana = gabwin(g,a_ana,kv.M);
g_syn = gabwin({'dual',{flags.wintype,gl,'inf'}},a_syn,kv.M);

%Initial coefficients
c = comp_sepdgtreal(f,g_ana,a_ana,kv.M,1);
s = abs(c);
phase = angle(c);

if flags.do_rtpghi || flags.do_pv || flags.do_phaselocked
    
    % Compute phase derivatives using either of the supported methods, see 
    % GABPHASEDERIV for more information. 
    if flags.do_dgt
        
        phased = gabphasederiv({'t','f'},'dgt',f,g_ana,a_ana,kv.M,'timeinv');
        
        % Sanitize and scale phase derivatives for further processing
        phased = ...
            cellfun(@(pEl) postpad(postpad(pEl,M2),newN,0,2),phased,'UniformOutput',0);
        [tgrad,fgrad] = deal(phased{:});
        tgrad = tgrad*2*pi/L;
        fgrad = fgrad*2*pi/L;
        
    elseif flags.do_phase
        
        % GABPHASEDERIV does not work natively with DGTREAL format 
        % coefficients. Therefore we compute derivatives directly.
        TimeInd = (0:(N-1))*a_ana;
        FreqInd = (0:(M2-1))/kv.M;
        phl = FreqInd'*TimeInd;
        phasefreqinv = phase-2*pi.*phl;
        tgrad = pderivunwrap(phasefreqinv,2,flags.diffmethod);
        tgrad = bsxfun(@plus,tgrad/a_ana,(0:M2-1)'*2*pi/kv.M);
        
        fgrad = pderivunwrap(phase,1,flags.diffmethod)/b_ana;
        fgrad(1,:) = fgrad(2,:);
        fgrad(end,:) = fgrad(end-1,:);
        
    elseif flags.do_abs
        
        % GABPHASEDERIV does not work natively with DGTREAL format 
        % coefficients. Therefore we replicate negative frequencies.
        s_full = [s;flipud(s(2:M2-1,:))];
        % The abs method assumes Gaussian window. Hence we pretend that the
        % analysis window is Gaussian. 
        gamma = pghi_findgamma(g,a_ana,kv.M,L);
        ggauss = {'gauss',gamma/L};
        phased = gabphasederiv({'t','f'},'abs',s_full,ggauss,a_ana,2,'timeinv');
        
        % Sanitize and scale phase derivatives for further processing
        phased = ...
            cellfun(@(pEl) postpad(postpad(pEl,M2),newN,0,2),phased,'UniformOutput',0);
        [tgrad,fgrad] = deal(phased{:});
        tgrad = tgrad*2*pi/L;
        fgrad = fgrad*2*pi/L;
    end
    
    
    % Process the DGT coefficients to realize the time stretching / 
    % compression.
    if flags.do_rtpghi % Phase Vocoder Done Right
        
        % RTPGHI accounts for scaling the phase derivatives according to
        % the synthesis parameters
        chat = rtpghi(s,{tgrad, fgrad},a_syn,kv.M,'timeinv');
        
    elseif flags.do_pv % Classical Phase Vocoder
        
        % Classical vocoder only considers time-direction derivative
        tgradint = tgrad*a_syn;
        newphase = cumtrapz([tgradint(:,end),tgradint,tgradint(:,1)],2);
        chat = s.*exp(1i*newphase(:,2:end-1));
        
    else
        error('%s: %s is not implemented yet.',upper(mfilename),flags.method);
    end

chat = postpad(chat,newN,0,2);
fout = postpad(comp_isepdgtreal(chat,g_syn,newL,a_syn,kv.M,1),exactnewL);

if flags.do_shift && stretch <= 1 % If pitch-shifting downwards, resample at the end
    fout = dctresample(fout,floor(exactnewL/stretch));  
end

end
info.a_syn = a_syn;
info.a_ana = a_ana;
info.M = kv.M;
info.c_ana = c;
info.c_syn = chat;

end 

function fd=pderivunwrap(f,dim,method)
%PDERIVUNWRAP Periodic derivative with unwrapping
%
%  `pderivunwrap(f,dim,wrapconst)` performs derivative of *f* along
%  dimmension *dim* using second order centered difference approximation
%  including unwrapping by *unwrapconst*
%
%  This effectivelly is just
%  (circshift(f,-1)-circshift(f,1))/2

if nargin<2 || isempty(dim)
    dim = 1;
end

unwrapconst = 2*pi;

shiftParam = {[1,0],[-1,0]};
if dim == 2
    shiftParam = {[0,1],[0,-1]};
end

switch method
    case 'centered'
        % Forward approximation
        fd_1 = f-circshift(f,shiftParam{1});
        fd_1 = fd_1 - unwrapconst*round(fd_1/(unwrapconst));
        % Backward approximation
        fd_2 = circshift(f,shiftParam{2})-f;
        fd_2 = fd_2 - unwrapconst*round(fd_2/(unwrapconst));
        % Average
        fd = (fd_1+fd_2)/2;
    case 'forward'
        fd = f - circshift(f,shiftParam{1});
        fd = fd - unwrapconst*round(fd/(unwrapconst));
    case 'backward'
        fd = circshift(f,shiftParam{2})-f;
        fd = fd - unwrapconst*round(fd/(unwrapconst));          
end


end 