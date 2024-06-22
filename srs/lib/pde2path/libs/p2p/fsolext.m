function [u,res,iter,Gu,Glam]=fsolext(p,u)
% FSOLEXT: replaces nloopext for using fsol.
%
%  [u,res,iter,Gu,Glam]=fsolext(p,u)
% Important settings: p.fsol.opt, p.fsol.tol, p.fsol.imax
%
% See also nloopext, fsol, fsolextf, stanparam
global pfs;
pfs=p; ua=u2au(p,u); 
opt=optimset(p.fsol.opt,'TolFun',p.fsol.tol,'MaxIter',p.fsol.imax);
[x,r,exitflag,output] = fsolve(@fsolextf,ua,opt);
u=au2u(p,x); res=norm(r,p.sw.norm); 
iter=output.iterations; [Gu,Glam]=getder(p,u,r);
end
