function Gu=getGuhpco(p,u,r)
% getGufoco: HP continuation extension of Gu.
p=hpreduce(p); % set regular case sizes to compute G and phi
nu=p.nu; M=p.mat.M(1:nu,1:nu); del=p.nc.del; ilam=p.nc.ilam; 
u1=[u(1:p.nu); u(3*p.nu+2:3*p.nu+1+p.naux)]; nu1=length(u1);  % org variables (no phi)
r1=r(1:nu); phr=u(nu+1:2*nu); phi=u(2*nu+1:3*nu); om=u(3*nu+1); % eigenvector and om 
Gv=getGupde(p,u1,r1); Gvr=Gv*phr; Gvi=Gv*phi;  % pde linearization and reference Gv*phi
if p.sw.spjac==1; Gvvph=p.fuha.spjac(p,u); % user def. 
else; Gvvph=sparse(2*nu,nu); % home-made expensive way
 for j=1:nu  
  up=u1+del*ej(j,nu1); rp=pderesi(p,up); 
  Gvp=getGupde(p,up,rp); % perturbed pde-part lin.
  rr=(Gvp*phr-Gvr)/del; ri=(Gvp*phi-Gvi)/del;   
  Gvvph=Gvvph+sparse(1:2*nu,j,[rr;ri],2*nu,nu); 
 end  
end
% first perturb by old primary parameter (other colums (nq>0) later) 
up=u1; up(nu+ilam(1))=up(nu+ilam(1))+del; 
rp=pderesi(p,up); Gv1=getGupde(p,up,rp); Gw=(rp-r1)/del;  
H2w=(Gv1*phr-Gvr)/del; H3w=(Gv1*phi-Gvi)/del; 
if(p.nc.nq==0)    
  Gu=[[ Gv                0*Gv   0*Gv   0*Gw   Gw ]; 
      [Gvvph(1:nu,:)      Gv     om*M   M*phi  H2w ]; 
      [Gvvph(nu+1:2*nu,:) -om*M  Gv    -M*phr  H3w ]; 
      [0*phr'             p.c  0*p.c   0     0  ];
      [0*phr'            0*p.c  p.c    0     0  ]]; 
  %figure(6); clf; spy(Gu), pause 
 % return; 
end
end