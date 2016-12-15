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
%         maxit   : Maximum number of iterations.
%         tol     : relative tolerance
%   Output parameters:
%         c       : Coefficients with the reconstructed phase
%         f       : Reconstructed signal.
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%
%   `decolbfgs(s,g,a,M)` attempts to find coefficients *c* from
%   their abs. value obtained as::
%     
%     c = dgtreal(f,g,a,M,'timeinv');
%     s = abs(c);
%
%   using Decorsiere's unconstrained optimization approach which, in turn,
%   exploits lBFGS (limited-memory Broyden–Fletcher–Goldfarb–Shanno)
%   algorithm. 
%
%   `[c,f,relres,iter] = decolbfgs(...)` additionally returns reconstructed
%   signal *f*, a vector of relative residuals *relres* and a number of
%   iterations done *iter*.
%
%   Initial phase estimate
%   ----------------------
%
%   'input'       Use phase of the input *s*. This is the default.
%
%   'zero'        Use zero phase.
%
%   'rand'        Use randomly generated phase.  
%
%   Additional parameters
%   ---------------------
%
%   'tol',tol     Tolerance for the lBFGS implementation. 
%                 The default value is 1e-6.
%
%   'maxit',maxit Maximum number of iterations of lBFGS.
%                 The default value is 100.
%
%   'p',p         p-parameter for the compressed version of the objective
%                 function. The default value is 2/3.
%   
%   'Ls',Ls       The reconstructed signal *f* will be cropped to length *Ls*.
%
%   'print'       Print info for each iteration. This is disabled by
%                 default.
%   
%   Note that this is just a wrapper around |frsynabs| from LTFAT.
%
%   Also note |minFunc| is required in order to run this function. It can be 
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


