function C=n2triamat(p,t)
% N2TRIAMAT: compute transf.-matrix from nodal values to triangle centers values  
%
%  C=n2triamat(p,t)
np=size(p,2); nt=size(t,2);
A=sparse(ones(3,1)*(1:nt),t(1:3,:),1,nt,np);
B=sparse(1:nt,1:nt,1./sum(A.'),nt,nt);
C=B*A; 
