function [Gua,Gun]=jaccheck(p)
% JACCHECK: compare numerical and assembled Jacobian (to check coding) 
% timing, and mesh quality.
r=resi(p,p.u); p.sw.jac=1;
tic; Gua=getGu(p,p.u,r); t2=toc; 
figure(p.plot.spfig2); spy(Gua); title('assembled G_{u}'); 
p.sw.jac=0;tic;Gun=getGu(p,p.u,r); t1=toc; 
figure(p.plot.spfig); spy(Gun); title('numerical G_{u}'); 
fprintf('time for numjac=%g, time for assembling=%g\n',t1,t2);
m1=full(max(max(abs(Gun)))); m2=full(max(max(abs(Gun-Gua)))); 
m3=full(max(sum(abs(Gun)))); m4=full(max(sum(abs(Gun-Gua))));
fprintf('max(Gun)=%g, max(Gun-Gua)=%g, infinity-norm(Gun)=%g, relerr=%g\n',...
    m1,m2,m3,m2/m3);
