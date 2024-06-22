function f=nodalf_schk(p,u) % for Schnakenberg 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end); % theta, gamma, d, a, b
x = getpte(p).';
theta = par(1);
f1 = p.myParams.betaf(x, theta) - u2.^2.*u1; 
f2 = u2.^2.*u1 - par(4) *u2 + p.myParams.eta(x,theta) ; 
f=[f1; f2];