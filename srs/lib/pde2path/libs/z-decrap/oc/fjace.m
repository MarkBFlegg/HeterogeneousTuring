function J=fjace(t,u) 
% fjace: Jac of rhs for OC BVP (in t) solver (iscarc)
% u=nodal vector(1:neq*np) at time t, 
% based on the discretization and p.jac contained in pdepath structure s1 
global s1 par; n=s1.nu; u=[u(1:n); par];  
if s1.sw.jac==0; r=resi(s1,u); else r=0; end; J=-getGupde(s1,u,r); 
J=[[J, zeros(n,1)]; zeros(1,n+1)]; 



