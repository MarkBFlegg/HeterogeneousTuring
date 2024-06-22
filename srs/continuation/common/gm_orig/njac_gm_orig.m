function [f1u,f1v,f2u,f2v]=njac_gm_orig(p,u) % Jacobian for Schnakenberg
%u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
par=u(p.nu+1:end); u1=u(1:p.nu/2); u2=u(p.nu/2+1:p.nu);

a_heterogeneity = p.myParams.a_heterogeneity;

% par=[theta gamma d a b];
theta = par(1);
gamma = par(2);
a = par(4);
b = par(5);

x = getpte(p)';

f1u = gamma * ( ...
    - b ...
    + 2 * u1 ./ u2 ...
);
f1v = gamma * ( ...
    - u1.^2 ./ (u2.^2) ...
);
f2u = gamma * ( ...
    2 * u1 ...
);
f2v = gamma * ( ...
    - 1 ...
) + 0.0*x;