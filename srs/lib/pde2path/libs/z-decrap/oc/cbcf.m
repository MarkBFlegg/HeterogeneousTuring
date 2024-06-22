function bc=cbcf(ya,yb) 
% cbcf: BC for OC with iscnat; require ya=u0 and psi*(yb-u1)=0 
global u0 u1 Psi; 
n=length(ya); bc=[ya(1:n/2)-u0(1:n/2); Psi*(yb-u1)]; 
