function [c,newphase,tgrad,fgrad]=rtpghifbwl(s,a,fc,tfr,varargin)
%RTPGHIFB Real-Time Phase Gradient Integration
%   Usage:  c=rtpghifb(s,gamma,a,fc);
%           [c,newphase,tgrad,fgrad] = rtpghifb(...);
%
%   Input parameters:
%         s        : Initial coefficients.
%         a        : Hop factor.
%         fc       : Vector of center frequencies.
%         tfr      : Time frequency ratio.
%   Output parameters:
%         c        : Coefficients with the constructed phase.
%         newphase : Just the (unwrapped) phase.
%         tgrad    : Relative time phase derivative.
%         fgrad    : Relative frequency phase derivative.
% 
%   `rtpghifb(s,gamma,a,fc)` creates filterbank coefficients from their
%   absolute values *s* using the Real-Time Phase Gradient Heap Integration
%   algorithm. *s* must have been obtained as::
%
%       c = ufilterbank(f,g,a,M,'timeinv');
%       s = abs(c);
%
%   and the algorithm attempts to recover *c*. Parameter *gamma* is window 
%   *g* specific and it can be computed using |pghi_findgamma|.
%
%   This function works entirely similar to |pghi| except it is using
%   the real-time version of the algorithm. Please see help of |pghi| 
%   (resp. |constructphasereal| from LTFAT) for more details.
%
%   Algorithm version:
%
%       'normal'    1 frame delay version of the algorithm RTPGHI(1)
%
%       'causal'    No delay version of the algorithm RTPGHI(0)
%
%   See also: dgtreal, idgtreal, pghi
%
%   References: 
%

%   AUTHORS: Clara Hollomey, Nicki Holighaus, Zdenek Prusa
%
thismfilename = upper(mfilename);
complainif_notposint(a,'a',thismfilename);
%complainif_notposint(fc,'fc',thismfilename);

% if ~( isscalar(gamma) || ...
%     ( iscell(gamma) && numel(gamma) == 2 ...
%     && all(cellfun(@(gEl) isequal(size(gEl),size(s)),gamma)) ...
%     && all(cellfun(@isreal,gamma))) )
%        error(['%s: gamma must be either a scalar or a 2 element cell array',...
%            ' containing phase derivatives such that {tgrad,fgrad}.'],...
%            upper(mfilename));
% end

definput.keyvals.tol=1e-6;
definput.flags.phase={'timeinv','freqinv'};
definput.flags.variant={'normal','causal'};
[flags,kv]=ltfatarghelper({},definput,varargin);
tol = kv.tol;
[M,N,W] = size(s);

if W>1
    error('%s: *s* must not be 3 dimensional.',thismfilename);
end

if numel(fc) ~= M
    error('%s: Mismatch between *fc* and the size of *s*.',thismfilename);
end

abss = abs(s);


[tgrad,fgrad] = ...
      comp_ufilterbankphasegradfrommag(...
      abss,N(1),a,M,tfr,fc,1,0);%flags.do_real...1, flags.do_gabor...0


logs=log(s + eps);
newphase = zeros(M,N);
c = zeros(M,N);

%scale fgrad and tgrad to match the normalized fc
fgrad = -pi*fgrad;
tgrad = a(1)*(tgrad + repmat(fc', 1,size(fgrad,2))) * pi;%scale to samples

for n=1:N
    idx = mod( n-1-1:n-1, N ) + 1;
    nprev = mod( n-2, N ) + 1;
    newphase(:,n) = comp_rtpghifbupdate(logs(:,idx),fc, tgrad(:,idx),fgrad(:,n),newphase(:,nprev),tol,(numel(fc)-1)*2);
    %newphase(:,n) = comp_rtpghifbupdate(abss(:,idx),info.fc,tgrad2(:, idx), fgrad2(:,n), newphase(:,nprev),tol,(numel(fc)-1)*2);
    
    c(:,n)=abss(:,n).*exp(1i*newphase(:,n));
end