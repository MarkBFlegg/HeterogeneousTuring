function r=sG_gm(p,u) % Schnakenberg 
f=nodalf_gm(p,u); par=u(p.nu+1:end);
Kscale=kron([[1,0];...
    [0,par(3)]],speye(size(p.mat.K,1)));
K = kron(speye(2),p.mat.K);
r=Kscale*(K*u(1:p.nu))-p.mat.M*f;