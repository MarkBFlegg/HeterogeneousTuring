function djca=disjca(p,sol,rho)
% disjca: discount jca (use for CSS)
par=p.u(p.nu+1:end); tl=length(sol.x); t=sol.x(tl); 
djca=jca(p,[sol.y(:,tl);par])*exp(-rho*t)/rho; 

