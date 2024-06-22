function [p] = interpolateFFT(p,nx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
oldnx = p.np;
par = p.u(p.nu+1:end);
lx = p.vol / 2;
p.pdeo=stanpdeo1D(lx,2*lx/nx);
p.np=p.pdeo.grid.nPoints;
p.nu=p.np*p.nc.neq; p.nc.neig=200;

nchemicals = p.nc.neq;
oldu = p.u;
p.u = zeros(nchemicals*p.np + numel(par),1);
for i=1:nchemicals
    p.u(nx * (i-1)+1:nx*i) = interpfft(oldu(oldnx* (i-1)+1:oldnx*i), nx);
end

p=setfemops(p);
end

