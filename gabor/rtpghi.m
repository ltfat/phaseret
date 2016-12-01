function [c,newphase,tgrad,fgrad]=rtpghi(s,gamma,a,M,varargin)
%RTPGHI Real-Time Phase Gradient Integration
%   Usage:  c=rtpghi(s,gamma,a,M);
%           c=rtpghi(c,gamma,a,M,mask);
%           c=rtpghi(s,gamma,a,M,mask,usephase);
%           [c,newphase,tgrad,fgrad] = rtpghi(...);
%
%   Input parameters:
%         s        : Initial coefficients.
%         gamma    : Window width factor.
%         a        : Hop factor.
%         M        : Number of channels.
%         tol      : Relative tolerance.
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
%       c = dgtreal(f,g,a,M);
%       s = abs(c);
%
%   and the algorithm attempts to recover *c*. Parameter *gamma* is window 
%   *g* specific and it can be computed using:
%
%   .. gamma = Cg*gl^2
%
%   .. math:: \gamma = C_g \mathit{gl}^2
%
%   where *gl* is the window length and *Cg* is a window specific constant. 
%   Both *Cg* and *gamma* can be obtained by calling |findwindowconstant|.
%
%   This function works entirely simiral to |pghi| except it is using
%   the real-time version of the algorithm. Please see help of |pghi| 
%   (resp. |constructphasereal| from LTFAT) for more details.
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

definput.keyvals.tol=1e-6;
definput.flags.phase={'timeinv','freqinv'};
definput.flags.variant={'normal','causal'};
[flags,kv]=ltfatarghelper({'mask','usephase'},definput,varargin);
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

[tgrad, fgrad, logs] = comp_pghiphasegrad( abss, gamma, a, M, flags.do_timeinv, flags.do_causal);

newphase = zeros(M2,N);
c = zeros(M2,N);

for n=1:N
    idx = mod( n-1-1:n-1, N ) + 1;
    nprev = mod( n-2, N ) + 1;
    newphase(:,n) = comp_rtpghiupdate(logs(:,idx),tgrad(:,idx),fgrad(:,n),newphase(:,nprev),tol,M);
    c(:,n)=abss(:,n).*exp(1i*newphase(:,n));
end
