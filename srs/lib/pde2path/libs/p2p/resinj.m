function r=resinj(t,au)
% RESINJ: dummy interface for residual via numerical jacobian numjac
%
%  r=resinj(t,au)
%
% See also numjac, pderesi, getGu
global pj 
p=pj;
ua=[au;p.u(p.nu+1:end)]; 
r=pderesi(p,ua); 
end
