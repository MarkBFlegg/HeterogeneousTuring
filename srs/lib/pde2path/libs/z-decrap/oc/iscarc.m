function [alv,vv,usec,esol,tlv,tv,uv]=iscarc(esol,usec,opt,fn)
% iscarc: Initial State Continuation ARClength (parametrization) 
% in : esol=starting point; usec=secant of last 2 two points 
%      opt=options for MTOM and add.options 
%      fn=filenames of CSS for startup, 
% out: alv,vv=alpha- and value-vector of successful solves 
%      usec,esol=last secant,sol (including alpha)    
%      tlv=number of grid points at each step
%      tv(i,tlv(i))=grid at step i
%      uv(i,n+1,tlv(i))=sol at step i  
global sig s1 xi um1 um2; 
if opt.start==1;  % 2 steps with iscnat to have the secant at alvin(end) 
  opt.msw=1; fprintf('iscarc, calling iscnat\n'); 
  [alv,vv,sol,ydat]=iscnat(opt.alvin,[],[],opt,fn); 
  n=s1.nu; al00=alv(end-1); al0=alv(end); % first 2 alphas 
  m=length(sol.x); esol=sol;  % prepare extended (by alpha) BVP-sol
  um2=[ydat.um2; al00*ones(1,m)]; um1=[ydat.um1; al0*ones(1,m)];
  delal=um1(n+1,1)-um2(n+1,1); xi=1/n; % last stepsize in alpha 
  usec=(um1-um2)/delal; sn=bxinorm(usec,xi); usec=usec/sn; % extended secant 
end
n=s1.nu; Me=[[s1.mat.M, zeros(n,1)]; [zeros(1,n) 1]]; % extend the mass matrix  
%opt=tomset(opt,'FJacobian',@fjace); opt=tomset(opt,'BCJacobian',@cbcjace);
opt.FJacobian=@fjace; opt.BCJacobian=@cbcjace;
alv=[]; vv=[]; tlv=[]; opt.M=Me; 
if opt.retsw==1 % prepare arrays for path returns; careful, can be LARGE!
  tv=zeros(opt.nsteps,opt.Nmax+1); uv=zeros(opt.nsteps,n+1,opt.Nmax+1); 
else tv=[]; uv=[]; 
end
for i=1:opt.nsteps
 ok=0; fprintf('i=%i, sig=%g\n ',i,sig); rho=s1.u(n+opt.rhoi); 
 while (sig>opt.sigmin) && (ok==0); 
  oldx=esol.x; esol.y=um1+sig*usec;  % predictor 
  [sol1,infv]=mtom(@mrhse,@cbcfe,esol,opt); info=sol1.err; fprintf('flag=%i, ',info); 
  if info==0; % sol. found
    ok=1; esol=sol1;tl=length(esol.x); al=esol.y(n+1,1); delal=al-um1(n+1,1); 
    alv=[alv al]; tlv=[tlv tl]; um2=um1; um1=esol.y; % append vals, update um1, um2
    if opt.retsw==1; tv(i,1:tl)=esol.x; uv(i,1:n+1,1:tl)=esol.y; end 
    if tl~=size(um2,2); um2=interp1(oldx,um2',esol.x); % interpolate um2 to new mesh
           um2=um2'; end  
    usec=um1-um2; sn=bxinorm(usec,xi); usec=usec/sn; % new secant 
    fprintf('al=%g, delta-al=%g\n', al, delal);   
    rsol=esol; rsol.y=esol.y(1:n,:); al=esol.y(n+1,1); % reduced soln (no alpha) 
    jcaval=jcai(s1,rsol,rho)+disjca(s1,rsol,rho); vv=[vv jcaval]; 
    if (infv.itnl<2 && sig<opt.sigmax/1.1) sig=sig*1.1; end % increase sig
  else if sig>opt.sigmin; sig=sig/2; fprintf('reducing sig to %g\n', sig);  
      else return; end
  end
 end
end
