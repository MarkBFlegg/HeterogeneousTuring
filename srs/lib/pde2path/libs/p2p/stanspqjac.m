function Qvvph=stanspqjac(p,u) 
% stanspqjac: 2nd derivative of aux. eqns, usually zero, i.e., linear q-eqn!
Qvvph=sparse(p.nc.nq,p.nu);