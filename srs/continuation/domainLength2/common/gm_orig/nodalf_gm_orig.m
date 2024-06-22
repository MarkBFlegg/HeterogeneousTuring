function f=nodalf_gm_orig(p,u) % for Schnakenberg 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end);


x = getpte(p).';

theta = par(1);
a = par(4);
b = par(5);
gamma = par(2);

f1 = gamma * ( ...
    u1 .^ 2 ./ u2  - b .* u1  + a * p.myParams.a_heterogeneity(theta, x) ...
); 
f2 = gamma * ( ...
    u1 .^ 2 - u2 ...
); 
f=[f1; f2];