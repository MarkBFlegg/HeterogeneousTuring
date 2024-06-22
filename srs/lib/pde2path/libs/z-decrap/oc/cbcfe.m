function bc=cbcfe(ya,yb) 
% cbcfe: (extended) BC for OC with iscarc; 
global s0 s1 Psi um1 um2 sig xi; n=length(ya)-1; 
v=(um1(:,1)-um2(:,1))'; vn=bxinorm(v,xi); v=v/vn; % last secant 
dv=ya-um1(:,1); % vector to last point 
addbc=xi*v(1:n)*dv(1:n)+(1-xi)*v(n+1)*dv(n+1)-sig; % Pseudo arclength
al=ya(n+1);  % current cont-par-value 
z0=al*s0.u(1:n/2)+(1-al)*s1.u(1:n/2); % left BC
bc=[ya(1:n/2)-z0; ... % left BC 
    Psi*(yb(1:n)-s1.u(1:n)); addbc]; 
