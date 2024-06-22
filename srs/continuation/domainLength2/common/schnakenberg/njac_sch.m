function [f1u,f1v,f2u,f2v]=njac_sch(p,u) % Jacobian for Schnakenberg
%u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
par=u(p.nu+1:end); u1=u(1:p.nu/2); u2=u(p.nu/2+1:p.nu);
f1u=par(2)*(-1+2*u1.*u2);
f1v=par(2) * u1.^2; 
f2u=par(2)*(-2*u1.*u2); 
f2v=par(2)*(-u1.^2); 