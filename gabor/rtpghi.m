function [c,newphase,tgrad,fgrad]=rtpghi(s,gamma,a,M,varargin)
%RTPGHI Real-Time Phase Gradient Integration
%   Usage:  c=rtpghi(s,gamma,a,M);
%           [c,newphase,tgrad,fgrad] = rtpghi(...);
%
%   Input parameters:
%         s        : Initial coefficients.
%         gamma    : Window width factor.
%         a        : Hop factor.
%         M        : Number of channels.
%   Output parameters:
%         c        : Coefficients with the constructed phase.
%         newphase : Just the (unwrapped) phase.
%         tgrad    : Relative time phase derivative.
%         fgrad    : Relative frequency phase derivative.
% 
%   `rtpghi(s,gamma,a,M)` creates complex DGTREAL coefficients from their
%   absolute values *s* using the Real-Time Phase Gradient Heap Integration
%   algorithm. *s* must have been obtained as::
%
%       c = dgtreal(f,g,a,M,'timeinv');
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
%   References: ltfatnote040 ltfatnote043
%

%   AUTHORS: Zdenek Prusa
%
thismfilename = upper(mfilename);
complainif_notposint(a,'a',thismfilename);
complainif_notposint(M,'M',thismfilename);

if ~( isscalar(gamma) || ...
    ( iscell(gamma) && numel(gamma) == 2 ...
    && all(cellfun(@(gEl) isequal(size(gEl),size(s)),gamma)) ...
    && all(cellfun(@isreal,gamma))) )
       error(['%s: gamma must be either a scalar or a 2 element cell array',...
           ' containing phase derivatives such that {tgrad,fgrad}.'],...
           upper(mfilename));
end

definput.keyvals.tol=1e-6;
definput.flags.phase={'timeinv','freqinv'};
definput.flags.variant={'normal','causal'};
[flags,kv]=ltfatarghelper({},definput,varargin);
tol = kv.tol;
[M2,N,W] = size(s);

if W>1
    error('%s: *s* must not be 3 dimensional.',thismfilename);
end

M2true = floor(M/2) + 1;

if M2true ~= M2
    error('%s: Mismatch between *M* and the size of *s*.',thismfilename);
end

abss = abs(s);

if iscell(gamma)
    L = N*a;
    b = L/M;
    tgrad = gamma{1}*a;
    fgrad = gamma{2}*b;
    logs=log(s + eps);
else
    [tgrad, fgrad, logs] = comp_pghiphasegrad( abss, gamma, a, M, flags.do_timeinv, flags.do_causal);
end

newphase = zeros(M2,N);
c = zeros(M2,N);

for n=1:N
    idx = mod( n-1-1:n-1, N ) + 1;
    nprev = mod( n-2, N ) + 1;
    newphase(:,n) = comp_rtpghiupdate(logs(:,idx),tgrad(:,idx),fgrad(:,n),newphase(:,nprev),tol,M);
    c(:,n)=abss(:,n).*exp(1i*newphase(:,n));
end
