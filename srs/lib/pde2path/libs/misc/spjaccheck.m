function [Ja,Jn]=spjaccheck(p)
% SPJACCHECK: compare numerical and assembled pa_u(G_u phi) (for foldcont) 
%
%   [Ja,Jn]=spjaccheck(p)
%
% "Ja" : assembled jacobian via spjac 
% "Jn" : numerical jacobian
%
% See also getGu, numjac, errcheck
r=resi(p,p.u); p=spreduce(p); u=p.u; nu=p.nu; 
u1=[u(1:nu);u(2*nu+1:2*nu+p.naux)]; % non-sp variables (no phi)
nu1=length(u1); % this is nu+p.naux (with nu halved)
r1=[r(1:nu);r(2*nu+1:2*nu+p.nc.nq)]; % normal variables
ph=u(nu+1:2*nu); % EVec 
tic; Ja=p.fuha.spjac(p,u); t1=toc; % fast way 
figure(p.plot.spfig); spy(Ja); title('user spjac'); 
fprintf('time for assembling=%g\n',t1);
% expensive way
Gv=getGupde(p,u1,r1); Gvph=Gv*ph; % pde linearization and reference product
Gvvph=sparse(nu,nu); tic; 
for j=1:nu  
  Gv1=getGupde(p,u1+p.nc.del*ej(j,nu1),r1); % perturbed pde-part lin.
  r2=(Gv1*ph-Gvph)/p.nc.del;   % add column with finite diff. approx.:
  Gvvph=Gvvph+sparse(1:nu,j,r2,nu,nu); 
  %figure(1); clf; spy(Gvvph),figure(2); clf; spy(r2), pause
end  
t2=toc; Jn=Gvvph; 
figure(p.plot.ifig); clf; spy(Jn); title('numerical spjac'); 
fprintf('time for FD appr=%g\n',t2);
m1=full(max(max(abs(Jn)))); m2=full(max(max(abs(Jn-Ja)))); 
m3=full(max(sum(abs(Jn)))); m4=full(max(sum(abs(Jn-Ja))));
fprintf('max(Jn)=%g, max(Jn-Ja)=%g, infinity-norm(Jn)=%g, relerr=%g\n',...
    m1,m2,m3,m4/m3);
