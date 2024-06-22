function f=hosrhs(t,u,p,par,T) 
% hosrhs: pad u with par, then call resi
%11, t,T,pause 
u=[u(1:p.nu); par]; p.t=t;  p.T=T; % to pass t,T to rhs in resi
f=resi(p,u); f=-T*f; 


