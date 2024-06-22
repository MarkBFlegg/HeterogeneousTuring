function p=hpcontini(varargin)
% HPCONTINI: initialization for Hopf-point continuation 
%  p=hpcontini(dir,fname,inewpar)
%  p=hpcontini(dir,fname,inewpar,newdir)    - set problem name to newdir
%  p=hpcontini(dir,fname,inewpar,newdir,ds) - also set stepsize to ds
%  p=hpcontini(dir,fname,inewpar,newdir,ds,aux) - add.arguments in aux
%
%  same with p, instead of dir/fname 
nnsw=0; resetsw=1; ds=0; 
if ischar(varargin{1}) % old syntax
   dir=varargin{1}; fname=varargin{2}; inewpar=varargin{3}; 
   if nargin>3; nname=varargin{4}; nnsw=1; end 
   if nargin>4; ds=varargin{5}; end 
   if nargin>5; aux=varargin{6}; end 
   p=loadp(dir,fname);    %generates also p.mat.drop and p.mat.fill
else p=varargin{1}; inewpar=varargin{2}; 
   if nargin>2; nname=varargin{3}; if ~isempty(nname); nnsw=1; end; end 
   if nargin>3; ds=varargin{4}; end 
   if nargin>4; aux=varargin{5}; end 
end 
if (p.sol.ptype~=3 && p.sol.ptype~=1)
    fprintf('\nWarning: not a HP. Set p.sw.spcont to 3 (HP) by hand;\n'); 
    fprintf('the Newton loop may converge to the desired point\n'); 
end
%if p.nc.nq>0; fprintf('HP-cont with nq>0 not yet implemented; add constraints to standard G, see cahn-hilliard demo 2.\n'); 
%   return; end 
fprintf('Initializing ... \n'); 
% Compute leading eigenvector
M=getM(p); M=[[M zeros(p.nu,p.nc.nq)]; [zeros(p.nc.nq,p.nu) eye(p.nc.nq)]]; 
Gu=getGu(p,p.u,resi(p,p.u));
%om0=p.nc.eigref(2); p.sol.muv, 
try eignr=aux.eignr; catch; eignr=length(p.nc.eigref); end 
om0=p.sol.muv(eignr,1); % om for finding Evec
evopts.disp=0; [phi1,mu1]=eigs(Gu,M,1,om0,evopts); 
fprintf('real(mu)=%g, imag(mu)=%g.\n',real(mu1),imag(mu1)); phi1=phi1/norm(phi1); % normalize eig-fct. 
phr=real(phi1); phi=imag(phi1); om=imag(mu1); p.c=2*phr'; 
% r1=0; r2=Gu*phr+om*M*phi; r3=Gu*phi-om*M*phr;  r4=p.c*phr+0*phi'*phi-1; r5=p.c*phi; 
% fprintf('r1=%g, r2=%g, r3=%g, r4=%g, r5=%g \n',norm(r1,'inf'),norm(r2,'inf'),norm(r3,'inf'),norm(r4,'inf'),norm(r5,'inf')); pause 
% Create new initial vector: format is u=[u; phi_r; phi_i; om; pars]
p.naux=length(p.u)-p.nu; % store number of auxiliary variables 
p.u=[p.u(1:p.nu); real(phi1(1:p.nu)); imag(phi1(1:p.nu)); imag(mu1); p.u(p.nu+1:p.nu+p.naux)];
% set new active variables including those for linear part (length=p.nc.nq+2+p.nc.nq)
p.nc.ilam=[inewpar; reshape(p.nc.ilam,length(p.nc.ilam),1)];
% set new equation numbers
p.nc.neq=3*p.nc.neq; p.nu=3*p.nu+1; p.nc.nq=2*p.nc.nq+1; 
p.sw.spcont=3; p.sol.restart=1;
p=tripledrop(p); % triple the size of p.mat.drop and p.mat.fill
fprintf('New active parameters and their values:\n'); printaux(p,1)
if ds~=0; p.sol.ds=ds; end % other data initialization
if resetsw; p=resetc(p); end % reset counters and branch
if nnsw; [p,ok]=setfn(p,nname); if ok~=1; q=p; return; end
else fprintf('warning: file name prefixes unchanged.\n'); end
p.sol.deta=0; p.sol.ineg=-1; p.sw.spcalc=0;
fname=[p.file.pname,sprintf('%i',p.file.count),'.mat'];save(fname,'p'); 
p.sw.spjac=1; p.fuha.spjac=@hpjac; 
p.sw.bifcheck=0; p.sw.spcalc=0; p.nc.del=1e-4; % set some switches 
