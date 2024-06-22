function jcaval=jcai(p,sol,rho) 
% jcai: jca along canonical path 
tv=sol.x; tl=length(tv); jcav=zeros(1,tl); par=p.u(p.nu+p.nc.nq+1:end); 
for i=1:tl
    jcav(i)=exp(-rho*tv(i))*jca(p,[sol.y(1:p.nu+p.nc.nq,i);par]);
end
jcaval=trapz(tv,jcav); 
