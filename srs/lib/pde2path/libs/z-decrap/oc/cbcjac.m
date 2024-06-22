function [ja,jb]=cbcjac(ya,yb)  
% cbcjac: Boundary jacobian for OC with iscnat; 
global Psi; 
n=length(ya); ja=sparse(n,n); jb=ja; 
for i=1:n/2; ja(i,i)=1; end
jb(n/2+1:n,:)=Psi(1:n/2,:); 
