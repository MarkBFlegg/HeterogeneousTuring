function f=mrhse(t,u) 
% mrhse: (extended) rhs for OC BVP (in t) solver
%
%  f=mrhse(t,u) 
% u=nodal vector(1:neq*np) at time t, in segment k 
% based on the discretization contained in pdepath structure p
global s1 par; u=[u(1:s1.nu);par]; f=-resi(s1,u); 
f=[f;0]; 


