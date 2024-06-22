function [Psi,muv,d,t1]=getPsi(p)
% getPsi: get matrix Psi for projection on E_u of target CSS
if ~isfield(p.mat,'M'); bc=p.fuha.bc(p,p.u); 
    [~,p.mat.M,~,~,p.mat.bcG,~,~]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,...
    0,1,zeros(p.nc.neq,1)); 
end
n=p.nu+p.nc.nq; p.nc.neig=n; % calc all Evals!
r=resi(p,p.u); Gu=getGu(p,p.u,r); % Jacobian at u1
ti1=tic; fprintf('getting Psi, '); 
[ineg,muv,EVs]=spcalcoc(Gu,p.mat.M,p);  % think M on lhs, hence solve gen. EVP 
ti1=toc(ti1); fprintf('done in %gs, ',ti1); 
d=ineg-n/2; musmin= muv(n/2+1); t1=10/abs(real(musmin)); % choose T=1/smallest stable EVal
fprintf('n/2=%i, d=%i, suggested T=%g\n',n/2, d, t1);
Psi=EVs(:,1:n/2)'; % the unstable directions (note u_t=-G(u)