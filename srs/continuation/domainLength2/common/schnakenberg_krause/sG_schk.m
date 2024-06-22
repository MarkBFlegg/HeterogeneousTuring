function r=sG_schk(p,u) % Schnakenberg 
f=nodalf_schk(p,u); par=u(p.nu+1:end); % theta, epsilon, d, alpha
Kscale=kron([[1*p.myParams.lScaling.^2,0];...
    [0,par(3)*p.myParams.lScaling.^2]],speye(size(p.mat.K,1)));
K = kron(speye(2)*par(2).^2,p.mat.K);
r=Kscale*(K*u(1:p.nu))-p.mat.M*f;