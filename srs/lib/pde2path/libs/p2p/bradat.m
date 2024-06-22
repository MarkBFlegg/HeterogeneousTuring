function bra=bradat(p)
% BRADAT: set default non-user-defined branch data 
%
%  bra=bradat(p)
%
% bra=[p.file.count; p.sol.ptype; p.sol.ineg; getlam(p); 
%      p.sol.err; L2-norm 1st component]
ineg=p.sol.ineg; 
if isfield(p,'hopf'); ho=p.hopf; 
    try ineg=ho.ind; end;   
end
ineg=max(ineg); u=p.u; 
M=getM(p); n1=floor(p.nu/p.nc.neq); l2=sqrt(u(1:n1)'*M(1:n1,1:n1)*u(1:n1)); 
bra=[p.file.count; p.sol.ptype; ineg'; getlam(p); p.sol.err; l2]; 
p.bra=bra;  