%%
% SPJACCHECK: compare numerical and assembled pa_u(G_u phi) (for spcont) 
%
%   [Ja,Jn]=spjaccheck(p)
%
% "Ja" : assembled jacobian via spjac 
% "Jn" : numerical jacobian
%
% See also getGu, numjac, errcheck
function [Ja,Jn]=spjaccheck(p, varargin)
global pj phij; 
njsw=1; if nargin>1; njsw=varargin{1}; end 
r=resi(p,p.u); 
[p,spcont,prim]=spreduce(p); % set regular case sizes and turn off spcont in this block
u=p.u; nu=p.nu; u1=[u(1:nu);u(2*nu+1:2*nu+p.naux)]; % non-sp variables (no phi)
%size(p.u), nu, size(u1), pause
nu1=length(u1); % this is nu+p.naux (with nu halved)
r1=[r(1:nu);r(2*nu+1:2*nu+p.nc.nq)]; % normal variables
ph=u(nu+1:2*nu); %om=u(2*nu+p.naux+1:end); % EVec 
tic; Ja=p.fuha.spjac(p,u); t1=toc; % fast way 
figure(p.plot.spfig); spy(Ja); title('user spjac'); 
fprintf('time for assembling=%g\n',t1);
% expensive way
Gv=getGupde(p,u1,r1); % get pde part linearization
Gvph=Gv*u(p.nu+1:2*p.nu); % take reference from residual
if njsw==1; tic; % use numjac 
  thresh=1e-10*p.nc.del*ones(p.nu,1); 
  phij=ph; pj=p;  M=getM(p); spmat=ones(p.nc.neq);
  np=nu/p.nc.neq; M=kron(spmat,M(1:np,1:np));  S=M^2>0; u=u(1:nu);
  [Gvvph,njfac,njG]=numjac('spnj',0,u(:),Gvph,thresh,[],0,S,[]); 
  t2=toc; 
else Gvvph=sparse(nu,nu); tic; % home-made 
  for j=1:nu  
   Gv1=getGupde(p,u1+p.nc.del*ej(j,nu1),r1); % perturbed pde-part lin.
   r2=(Gv1*ph-Gvph)/p.nc.del;   % add column with finite diff. approx.:
   Gvvph=Gvvph+sparse(1:nu,j,r2,nu,nu); % Gvvph(:,j-1)=(r1-r)/del; 
  end  
  t2=toc; 
end
Jn=Gvvph; 
figure(p.plot.ifig); clf; spy(Jn); title('numerical spjac'); 
fprintf('time for FD appr=%g\n',t2);
m1=full(max(max(abs(Jn)))); m2=full(max(max(abs(Jn-Ja)))); 
m3=full(max(sum(abs(Jn)))); m4=full(max(sum(abs(Jn-Ja))));
fprintf('max(Jn)=%g, max(Jn-Ja)=%g, infinity-norm(Jn)=%g, relerr=%g\n',...
    m1,m2,m3,m4/m3);
end 