function f=nodalf_gm(p,u) % for Schnakenberg 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end);


x = getpte(p).';

theta = par(1);
a = par(4);
c = par(5);
p_param = par(6);
gamma = par(2);

f1 = gamma * ( ...
    u1 .^ 2 ./ ( (1 + p_param * u1.^2) .* u2 ) - c * p.myParams.c_heterogeneity(theta, x) .* u1 ...
); 
f2 = gamma * ( ...
    u1 .^ 2 - a * u2 .* p.myParams.a_heterogeneity(theta, x) ...
); 
f=[f1; f2];