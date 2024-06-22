function stancssvalf(dir,pt)  
% stancssvalf: print CSS value and characteristics
p=loadp(dir,pt); fprintf([dir '/' pt ', lam=%g, '],getlam(p)); 
r=p.u(p.nu+1); n=p.np; [pt,tr]=getpte(p); 
c=p.fuha.con(p,p.u); u1=p.u(1:n); ca=triint(c,pt,tr)/p.vol; ua=triint(u1,pt,tr)/p.vol; 
fprintf('(<u1>,<k>,jca)=(%4.2f & %4.2f & %4.2f)\n', ua,ca,jca(p,p.u)/r);


