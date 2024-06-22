function [ja,jb]=cbcjace(ya,yb)  
% cbcjace: Boundary jacobian for OC with iscarc (extended) 
global s0 s1 Psi um1 um2 xi; n=length(ya)-1; 
v=(um1(:,1)-um2(:,1))'; vn=bxinorm(v,xi); v=v/vn; 
ja=sparse(n+1,n+1); for i=1:n/2; ja(i,i)=1; end % left BC 
ja(1:n/2, n+1)=-s0.u(1:n/2)+s1.u(1:n/2); 
ja(n+1,1:n)=xi*v(1:n); ja(n+1,n+1)=(1-xi)*v(n+1); 
jb=sparse(n+1,n+1);jb(n/2+1:n,1:n)=Psi(1:n/2,1:n); % right BC

