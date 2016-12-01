function [c,f,relres,iter]=decolbfgs(s,g,a,M,varargin)
%DECOLBFGS Decorsiere's phase reconstruction algorithm
%   Usage: c = decolbfgs(s,g,a,M)
%          c = decolbfgs(s,g,a,M,maxit)
%          c = decolbfgs(s,g,a,M,maxit,tol)
%          [c,f,relres,iter] = decolbfgs(...)
%
%   Input parameters:
%         s       : Initial coefficients.
%         g       : Analysis Gabor window
%         a       : Hop factor
%         M       : Number of channels
%         Ls      : Length of signal.
%   Output parameters:
%         f       : Signal.
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%         c       : Coefficients with the reconstructed phase
%
%   `decolbfgs(s,g,a,M)` attempts to find coefficients *c*
%   their abs. value is::
%
%     s = abs(c)
%
%   using Decorsiere's unconstrained optimization approach which in turn
%   exploits lBFGS (limited-memory Broyden–Fletcher–Goldfarb–Shanno)
%   algorithm. 
%
%   Note |minFunc| is required in order to run this function. It can be 
%   obtained from here https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
%
%   References: desomada15
%

% AUTHOR: Zdenek Prusa
%

definput.flags.phase={'timeinv','freqinv'};
definput.keyvals.p=2/3;
definput.keyvals.tol=1e-6;
definput.keyvals.maxit=100;
definput.keyvals.Ls=[];
definput.keyvals.printstep=10;
definput.flags.print={'quiet','print'};
definput.flags.startphase={'input','zero','rand'};
[flags,kv]=ltfatarghelper({'maxit','tol'},definput,varargin);
Ls = kv.Ls;

F = frame('dgtreal',g,a,M,flags.phase);

sframe = framenative2coef(F,s);

[f,relres,iter] = frsynabs(F,sframe,'bfgs','p',kv.p,'tol',kv.tol,'maxit',kv.maxit,...
    flags.startphase,flags.print,'printstep',kv.printstep);

c = framecoef2native(F,frana(F,f));

if isempty(Ls)
    Ls = size(s,2)*a;
end

f = postpad(f,Ls);


