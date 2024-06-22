function [p] = init_gm(opts)
%INIT_SCH Summary of this function goes here
%   Detailed explanation goes here
p=[];

p.myParams.a_heterogeneity = opts.par.a_heterogeneity;
p.myParams.c_heterogeneity = opts.par.c_heterogeneity;


theta = opts.par.theta;
gamma = opts.par.gamma;
d = opts.par.d;
a = opts.par.a;
c = opts.par.c;
s = opts.par.s;
p_param = a/c/s/s/s - 1/s/s


%try to start at patterned steady state.
par=[theta gamma d a c p_param];

p=gminit(p,opts.L,opts.Nx,par);
p.sw.bifcheck=2;
 
x = getpte(p).';

u = s + x .* 0.0;
v = s*s/a + x.*0.0;
p.u=[u;v;par']; % initial solution guess with parameters
p = setfn(p, opts.fname);

end

