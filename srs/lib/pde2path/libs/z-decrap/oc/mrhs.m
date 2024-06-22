function f=mrhs(t,u,k) 
% mrhs: rhs for OC BVP (in t) solver
%
%  f=mrhs(t,u,k) 
% u=nodal vector(1:neq*np) at time t, in segment k 
% based on the discretization contained in pdepath structure s1
global s1 par; u=[u;par]; f=-resi(s1,u); 
%n=s1.np; f(1:3)', f(n+1:n+3)', pause 


