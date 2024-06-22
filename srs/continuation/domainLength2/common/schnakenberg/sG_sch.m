function r=sG_sch(p,u) % Schnakenberg 
f=nodalf_sch(p,u); par=u(p.nu+1:end); % theta, gamma, d, a, b
K=kron([[1,0];[0,par(3)]],p.mat.K); 
r=K*u(1:p.nu)-p.mat.M*f;