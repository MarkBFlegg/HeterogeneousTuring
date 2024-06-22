close all; clear all;
opts.fname = 'test';

opts.par.theta = 0.0;
opts.par.d = 20;
opts.par.gamma = 50000;
opts.par.a = 0.9;
opts.par.c = 1;
opts.par.s = 0.708;

opts.par.a_heterogeneity = ...
    @(theta,x) (1 + theta *cos(pi .* (x+0.5)));
opts.par.c_heterogeneity = ...
    @(theta,x) 2 - opts.par.a_heterogeneity(theta,x);

opts.L = 0.5;
opts.Nx = 1e3;

p = init_gm(opts);

% Set continuation parameters

p.nc.ilam=1;
p.nc.dsmin = 1e-5;
p.sol.ds=1e-6;
p.nc.dsmax=1e-1;
p.nc.dlammin = 1e-5;
p.nc.dlammax = 1e-1;
p.nc.lammin=0;
p.nc.lammax=1;

p.sw.eigssol = 0;
p.sw.bifcheck = 1;
p.sw.spcalc = 1;
p.nc.tol = 1e-8;
p.fsol.tol = 1e-10;
p.fsol.fsol = 0;

p.fsol.disp = 0;

p = cont(p, 1000);
%p = cont(p, 2000)

%p = cont(p, 100)


