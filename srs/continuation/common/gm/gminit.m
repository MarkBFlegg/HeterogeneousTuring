function p=gminit(p,dom,nx,par,varargin)
p.fuha = {};
p=stanparam(p);
%screenlayout(p); 
p.nc.neq=2; p.sw.sfem=-1;
p.sw.spjac=0; % do not use analytical Jacobian for spectral point cont (fold cont)
% theta gamma d a c p
p.plot.auxdict={'\theta','\gamma','d','a','c','p','||u_1||_{2}','min(|u_1|)'};
p.fuha.sG=@sG_gm;
p.fuha.sGjac=@sGjac_gm;

p.fuha.lss = @lss;%@myIterativeLss;
p.fuha.blss = @lss;%@myIterativeLss;
p.fuha.innerlss = @lss;%@myIterativeLss;
%p.fuha.spjac=@spjac; 
switch length(dom)
  case 1; lx=dom; p.pdeo=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; 
  case 2; ny=round(dom(2)/dom(1)*nx); lx=dom(1); ly=dom(2); p.vol=4*lx*ly; 
      pde=stanpdeo2D(lx,ly,nx,ny,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
  case 3; lx=dom(1); ly=dom(2); lz=dom(3); h=2*lx/(nx-1);  p.vol=8*lx*ly*lz; 
      pde=stanpdeo3D(lx,ly,lz,h,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
      p.plot.EdgeColor='none'; 
end
p.np=p.pdeo.grid.nPoints; p.nu=p.np*p.nc.neq; p.nc.neig=20; p.nc.nsteps=5000; 
p=setfemops(p); p.sol.xi=1/p.nu; p.file.smod=10; p.sw.para=2; p.sw.foldcheck=1; 
p.sw.norm=2;

p.fsol.fsol = 0;

u=ones(p.np,1);
v=ones(p.np,1); 

p.u=[u;v;par']; % initial solution guess with parameters
[po,tr,ed]=getpte(p); p.mesh.bp=po; p.mesh.be=ed; p.mesh.bt=tr; p.nc.ngen=1; 
p.plot.bpcmp=length(par)+1; p.plot.axis='image'; p.plot.cm='hot'; 
p.nc.resfac=1e-3; p.pm.mst=4; 