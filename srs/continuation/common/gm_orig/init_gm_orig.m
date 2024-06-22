function [p] = init_gm_orig(opts)
%INIT_SCH Summary of this function goes here
%   Detailed explanation goes here
p=[];

p.myParams.a_heterogeneity = opts.par.a_heterogeneity;

theta = opts.par.theta;
gamma = opts.par.gamma;
d = opts.par.d;
a = opts.par.a;
b = opts.par.b;


%try to start at patterned steady state.
par=[theta gamma d a b];

p=gm_originit(p,opts.L,opts.Nx,par);
p.sw.bifcheck=2;
 
x = getpte(p).';

u = (1+a)/b + x * 0;
v = u.*u;
p.u=[u;v;par']; % initial solution guess with parameters
p = setfn(p, opts.fname);

end

