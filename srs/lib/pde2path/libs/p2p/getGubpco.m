function Gu=getGubpco(p,u,r)
% getGubpco: BP continuation extension of Gu. 
% The assumption is that initially this will be called after bpcontini, 
% so that the old prim. param. is at p.nc.ilam(2).
%
%  v=upde normal unknowns part,   w=normal active uaux part
%  ph=u(p.nu/2+1:p.nu) upde part of adjoint kernel eigenvector
%  om=u(p.nu+p.naux+1:length(u)) active non-primary uaux part of kernel eigenvector 
%  G=normal pde equations, Q=normal auxiliary equations
%  Gv=partial derivative G_v, others analogous
%  Gvvph=\pa_v(G_v'*ph)
ilam=p.nc.ilam; 
p=bpreduce(p); % set regular case sizes to compute G and Gu'*phi
u1=[u(1:p.nu); u(2*p.nu+1:2*p.nu+p.naux)]; % non-sp variables (no psi, mu)
nu1=length(u1); % this is p.nu+p.naux (with p.nu halved)
r1=[r(1:p.nu); r(2*p.nu+1:2*p.nu+p.nc.nq)]; % normal variables
ph=u(p.nu+1:2*p.nu); om=u(2*p.nu+p.naux+2:end); mu=u(2*p.nu+ilam(end)); 
Gv=getder(p,u1,r1); % get analytic pde part linearization Glam=wrt old lam
Glam=getGlam(p,u1,r1); 
Gvph=Gv'*ph; % take reference from residual  (ph=psi here!!!) 
if p.sw.spjac==1; Gvvph=p.fuha.spjac(p,u); % user def. 
else % expensive way
 Gv=getGupde(p,u1,r1); % get pde part linearization
 Gvph=Gv*u(p.nu+1:2*p.nu); % take reference from residual
 Gvvph=sparse(p.nu,p.nu); % home-made 
 for j=1:p.nu  
   Gv1=getGupde(p,u1+p.nc.del*ej(j,nu1),r1); % perturbed pde-part lin.
   r2=(Gv1-Gv)'*ph/p.nc.del;
   Gvvph=Gvvph+sparse(1:p.nu,j,r2,p.nu,p.nu); % add column with finite diff. approx.:
 end   
end
% finite difference for mixed pde-aux part Gvwph=\pa_w(G_v*ph) and Gwvom=\pa_w(G_a*om)
uper=u1+p.nc.del*ej(p.nu+p.nc.ilam(1),nu1); % perturb by old primary parameter
r2=pderesi(p,uper); Gv2=getGupde(p,uper,r2);
% part of mixed deriv. multiplied with ph (pde part of evec.)
Gvwph=(Gv2-Gv)'*ph/p.nc.del;  Gw=(r2-r(1:p.nu))/p.nc.del; 
if(p.nc.nq==0)  
  % assemble Gu from block matrices; 3rd row is the normalization
  M=p.mat.M(1:p.nu,1:p.nu); 
  Gu=[[ Gv   mu*M  Gw    M*ph];  % G+mu*psi
      [Gvvph   Gv' Gvwph  0*ph]; % G^T*psi
      [0*ph' 2*ph'   0    0]     % |psi|^2-1
      [Gvwph'  Glam'  0   0]];   % Glam'*psi 
 % figure(6); clf; spy(Gu), pause 
  return; 
end
% now case nq>0; exclude new primary parameter, and initially also old one:
 Gwvom=sparse(p.nu,p.nu); 
for k=1:p.nc.nq 
       Gv1=getGupde(p,u1+p.nc.del*ej(p.nu+p.nc.ilam(k+1),nu1),r1); 
       % part of mixed deriv. multiplied with ph (pde part of evec.)
       r2=(Gv1*ph-Gvph)/p.nc.del;
       Gvwph=Gvwph+sparse(1:p.nu,k+1,r2,p.nu,p.nc.nq+1);
       % part of mixed deriv. multiplied with components of om (aux part of evec.)
       Gwvom=Gwvom+u(2*p.nu+p.naux+k)*(Gv1-Gv)/p.nc.del;
end
% case of (org) nq>0; ****************** Needs revision !!!! 
Qvvph=sparse(p.nc.nq,p.nu); Qwvom=sparse(p.nc.nq,p.nu); Qvwph=sparse(p.nc.nq,p.nc.nq+1);
rq=r(2*p.nu+1:2*p.nu+p.nc.nq);       % aux. eqn. residual
Qvph=p.fuha.qfder(p,u1)*ph;          % aux. eqn. linearization residual
% this is in fact also r(2*p.nu+p.nc.nq+1:length(r)-1)', but an initially poor guess 
% might cause divergence
 if p.sw.qjac==1 % analytical q_u 
   Qv=p.fuha.qfder(p,u1); spqjac=0; if isfield(p.sw,'spqjac'); spqjac=p.sw.spqjac; end 
   if spqjac==1; Qvvph=p.fuha.spqjac(p,u);  % analytical pa_u(q_u*phi) 
   else
     for j=1:p.nu
       Qv1=p.fuha.qfder(p,u1+p.nc.del*ej(j,nu1)); % perturbed pde-part lin.
       r2=(Qv1*ph-Qvph)/p.nc.del;
       Qvvph=Qvvph+sparse(1:p.nc.nq,j,r2,p.nc.nq,p.nu);
     end
   end
   Qv1=p.fuha.qfder(p,u1+p.nc.del*ej(p.nu+p.nc.ilam(1),nu1)); 
   r2=(Qv1*ph-Qvph)/p.nc.del;
   Qvwph=Qvwph+sparse(1:p.nc.nq,1,r2,p.nc.nq,p.nc.nq+1); % parameter der. of q_u*phi
   for k=1:p.nc.nq
       Qv1=p.fuha.qfder(p,u1+p.nc.del*ej(p.nu+p.nc.ilam(k+1),nu1)); 
       r2=(Qv1*ph-Qvph)/p.nc.del;
       Qvwph=Qvwph+sparse(1:p.nc.nq,k+1,r2,p.nc.nq,p.nc.nq+1);
       Qwvom=Qwvom+u(2*p.nu+p.naux+k)*(Qv1-Qv)/p.nc.del;
   end
else % all second order numerically
   Qv=sparse(p.nc.nq,p.nu); 
   for j=1:p.nu
       r1=p.fuha.qf(p,u1+p.nc.del*ej(j,nu1));
       Qv=Qv+sparse(1:p.nc.nq,j,(r1-rq)/p.nc.del,p.nc.nq,p.nu);           
       for jj=1:p.nu  % second derivatives: SLOW!
          rpp=p.fuha.qf( p,u1+p.nc.del*(ej(j,nu1)+ej(jj,nu1)) ); 
          rmm=p.fuha.qf( p,u1-p.nc.del*(ej(j,nu1)+ej(jj,nu1)) ); 
          rpm=p.fuha.qf( p,u1+p.nc.del*(ej(j,nu1)-ej(jj,nu1)) );
          rmp=p.fuha.qf( p,u1-p.nc.del*(ej(j,nu1)-ej(jj,nu1)) );
          r2=u(p.nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
          Qvvph=Qvvph+sparse(1:p.nc.nq,j,r2,p.nc.nq,p.nu);
       end
       %   mixed derivatives:
       %   first with respect to in initial run the old primary parameter
       kk=p.nu+p.nc.ilam(1);
       rpp=p.fuha.qf( p,u1+p.nc.del*(ej(j,nu1)+ej(kk,nu1)) ); 
       rmm=p.fuha.qf( p,u1-p.nc.del*(ej(j,nu1)+ej(kk,nu1)) ); 
       rpm=p.fuha.qf( p,u1+p.nc.del*(ej(j,nu1)-ej(kk,nu1)) );
       rmp=p.fuha.qf( p,u1-p.nc.del*(ej(j,nu1)-ej(kk,nu1)) );
       r2=u(p.nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
       Qvwph=Qvwph+sparse(1:p.nc.nq,1,r2,p.nc.nq,p.nc.nq+1);
       %   now with respect to other old active aux. vars
       for k=2:p.nc.nq+1 
          kk=p.nu+p.nc.ilam(k);
          rpp=p.fuha.qf( p,u1+p.nc.del*(ej(j,nu1)+ej(kk,nu1)) ); 
          rmm=p.fuha.qf( p,u1-p.nc.del*(ej(j,nu1)+ej(kk,nu1)) ); 
          rpm=p.fuha.qf( p,u1+p.nc.del*(ej(j,nu1)-ej(kk,nu1)) );
          rmp=p.fuha.qf( p,u1-p.nc.del*(ej(j,nu1)-ej(kk,nu1)) );
          r2=u(2*p.nu+p.naux+k-1)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
          Qwvom=Qwvom+sparse(1:p.nc.nq,j,r2,p.nc.nq,p.nu);
          r2=u(p.nu+j)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
          Qvwph=Qvwph+sparse(1:p.nc.nq,k-1,r2,p.nc.nq,p.nc.nq+1);
       end
   end
end % qfder=1  
% Notation: 'lam' here means the (initially) old primary parameter
% In case p.nc.nq=0 above this was named 'Gw' while here Gw stands for 
% other aux. vars. derivatives.
Gwlam=sparse(p.nu,p.nc.nq+1); Gwwom=sparse(p.nu,p.nc.nq+1);
Qwlam=sparse(p.nc.nq,p.nc.nq+1); Qwwom=sparse(p.nc.nq,p.nc.nq+1);
for k=2:p.nc.nq+2
    k1=p.nu+p.nc.ilam(k-1);
       % Gw finite diff via pde:
       r1=pderesi(p,u1+p.nc.del*ej(k1,nu1));
       Gwlam=Gwlam+sparse(1:p.nu,k-1,(r1-r(1:p.nu))/p.nc.del,p.nu,p.nc.nq+1);           
       % Qw finite diff:
       r1=p.fuha.qf(p,u1+p.nc.del*ej(k1,nu1));
       Qwlam=Qwlam+sparse(1:p.nc.nq,k-1,(r1-rq)/p.nc.del,p.nc.nq,p.nc.nq+1);
       % second derivatives:
       for kk=3:p.nc.nq+2 % start at 2 to skip primary parameter
          kk1=p.nu+p.nc.ilam(kk-1);
          rpp=p.fuha.qf( p,u1+p.nc.del*(ej(k1,nu1)+ej(kk1,nu1)) );
          rmm=p.fuha.qf( p,u1-p.nc.del*(ej(k1,nu1)+ej(kk1,nu1)) );
          rpm=p.fuha.qf( p,u1+p.nc.del*(ej(k1,nu1)-ej(kk1,nu1)) );
          rmp=p.fuha.qf( p,u1-p.nc.del*(ej(k1,nu1)-ej(kk1,nu1)) );
          r2=u(2*p.nu+p.naux+kk-2)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
          Qwwom=Qwwom+sparse(1:p.nc.nq,k-1,r2,p.nc.nq,p.nc.nq+1);
          % now for Gww... 
          rpp=pderesi( p,u1+p.nc.del*(ej(k1,nu1)+ej(kk1,nu1)) ); 
          rmm=pderesi( p,u1-p.nc.del*(ej(k1,nu1)+ej(kk1,nu1)) ); 
          rpm=pderesi( p,u1+p.nc.del*(ej(k1,nu1)-ej(kk1,nu1)) );
          rmp=pderesi( p,u1-p.nc.del*(ej(k1,nu1)-ej(kk1,nu1)) );
          r2=u(2*p.nu+p.naux+kk-2)*(rpp-rmp-rpm+rmm)/(4*p.nc.del^2);
          Gwwom=Gwwom+sparse(1:p.nu,k-1,r2,p.nu,p.nc.nq+1);
       end
   end
   % remove spcont-parameter; initially the old primary parameter:
   Gw=Gwlam(:,2:end); 
   Qw=Qwlam(:,2:end); 

   % assemble Gu from block matrices; last row is the normalization
   zeroGv=sparse(p.nu,p.nu); zeroGw=sparse(p.nu,p.nc.nq); 
   zeroQw=sparse(p.nc.nq,p.nc.nq); zeroQv=sparse(p.nc.nq,p.nu);
   zerovt=sparse(1,p.nu);    zerowt=sparse(1,p.nc.nq+1);
   Gu=[[    Gv      zeroGv    Gwlam     zeroGw];
        [Gvvph+Gwvom   Gv   Gvwph+Gwwom    Gw  ];
        [    Qv      zeroQv    Qwlam     zeroQw];
        [Qvvph+Qwvom   Qv   Qvwph+Qwwom    Qw  ];
        [  zerovt     2*ph'    zerowt     2*om']];
end 
