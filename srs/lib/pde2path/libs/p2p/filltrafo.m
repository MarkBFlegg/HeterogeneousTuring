function M=filltrafo(p,M)
% FILLTRAFO: return M transformed via p.mat.fill to change M from Neumann b.c. to periodic b.c..
%  (unless p.mat.fill=1 in which case the matrix M is left unchanged).
%
%  M=filltrafo(p,M)
%
% See also getCylOp, getTorOp
n=size(M,1); mdim=round(n/p.np); %size(M), pause 
if size(p.mat.fill,1)>1 
    fM=p.mat.fill(1:mdim*p.np,1:mdim*p.nu/p.nc.neq); 
%M=p.mat.fill'*M*p.mat.fill; 
    M=fM'*M*fM; 
end 
%M=p.mat.drop*M*p.mat.fill; 