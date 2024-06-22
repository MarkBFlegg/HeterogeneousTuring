function [p] = init_sch(opts)
%INIT_SCH Summary of this function goes here
%   Detailed explanation goes here
p=[];

a = opts.par.a;
b = opts.par.b;
u1 = opts.par.a+opts.par.b;
u2 = b/((a+b)^2);

%try to start at patterned steady state.
par=[0,opts.par.gamma,opts.par.d, opts.par.a, ...
    opts.par.b]; % theta, gamma, d, a, b

p=schnakinit(p,opts.L,opts.Nx,par);
p.sw.bifcheck=2; 

x = getpte(p).';

u = u1 + (x .* 0);
v = u2 + (x .* 0);
p.u=[u;v;par']; % initial solution guess with parameters
p = setfn(p, opts.fname);

end

