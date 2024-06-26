function Gu=sGjac_schk(p,u)
par=u(p.nu+1:end); n=p.np;
[f1u,f1v,f2u,f2v]=njac_schk(p,u); % (nodal) Jacobian of 'nonlin' 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=kron([[1*p.myParams.lScaling.^2,0];[0,par(3)*p.myParams.lScaling.^2]],p.mat.K.*par(2).^2)-p.mat.M*Fu; 