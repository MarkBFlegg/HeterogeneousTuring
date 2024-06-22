function [p] = init_schk(opts)
%INIT_SCH Summary of this function goes here
%   Detailed explanation goes here
p=[];

alpha = opts.par.alpha;
betaf = opts.par.betaf;
p.myParams.betaf = betaf;
p.myParams.eta = opts.par.eta;
theta = opts.par.theta;

%try to start at patterned steady state.
par=[theta,opts.par.epsilon,opts.par.d, opts.par.alpha];

p=schnakinitk(p,opts.L,opts.Nx,par);
p.sw.bifcheck=2;

x = getpte(p).';

u = alpha * alpha * betaf(x,theta);
v = 1 ./ alpha + x.*0.0;
p.u=[u;v;par']; % initial solution guess with parameters
p = setfn(p, opts.fname);

end

