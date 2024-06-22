function p=setbmesh(p)
% SETBMESH: (re)set basemesh for meshadac to current mesh.
%
%  p=setbmesh(p) 
%
% See also meshadac, stanparam.
p.mesh.bp=p.mesh.p; p.mesh.be=p.mesh.e; p.mesh.bt=p.mesh.t; p.mesh.maxt=2*p.mesh.nt; 
