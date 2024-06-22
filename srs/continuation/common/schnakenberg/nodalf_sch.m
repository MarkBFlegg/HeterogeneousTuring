function f=nodalf_sch(p,u) % for Schnakenberg 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end); % theta, gamma, d, a, b
x = getpte(p).';
theta = par(1);
f1=par(2)*(par(4)-u1+u1.^2.*u2 + p.myParams.gu(x, theta)); 
f2=par(2) * (par(5)-u1.^2.*u2 + p.myParams.gv(x, theta)); 
f=[f1; f2]; 