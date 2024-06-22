function [p] = interpolateCub(p,nx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
oldnx = p.np;
oldx = getpte(p).';
par = p.u(p.nu+1:end);
lx = p.vol / 2;
p.pdeo=stanpdeo1D(lx,2*lx/nx);
p.np=p.pdeo.grid.nPoints;
p.nu=p.np*p.nc.neq; p.nc.neig=200;
newnx = numel(getpte(p));

nchemicals = p.nc.neq;
oldu = p.u;
p.u = zeros(nchemicals*p.np + numel(par),1);

olddx = oldx(2) - oldx(1);
oldx = [oldx(1) - olddx; oldx; oldx(end) + olddx];
for i=1:nchemicals
    interpu = oldu(oldnx* (i-1)+1:oldnx*i);
    interpu = [interpu(2); interpu; interpu(end-1)];
    p.u(newnx * (i-1)+1:newnx*i) = spline(oldx,...
        interpu, getpte(p));
end
p.u(p.nu+1:end) = par;

p=setfemops(p);
end

