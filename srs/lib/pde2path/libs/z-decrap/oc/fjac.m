function J=fjac(t,u) 
% fjac: Jac of rhs for OC BVP (in t) 
% u=nodal vector(1:neq*np) at time t, 
% based on the discretization and p.jac contained in pdepath structure s1
global s1 par; u=[u; par]; 
if s1.sw.jac==0; r=resi(s1,u); else r=0; end; J=-getGu(s1,u,r); 




