function [f1u,f1v,f2u,f2v]=njac_schk(p,u) % Jacobian for Schnakenberg
%u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
par=u(p.nu+1:end); u1=u(1:p.nu/2); u2=u(p.nu/2+1:p.nu);
f1u= (-u2.^2);
f1v= (-2*u1.*u2);
f2u= (u2.^2);
f2v= (2* u1.*u2 - par(4));