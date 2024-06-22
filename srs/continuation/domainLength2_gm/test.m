close all; clear all;
opts.fname = 'test';

opts.par.theta = 0.0;
opts.par.d = 20;
opts.par.a = 0.1;
opts.par.b = 1;

u_star = (opts.par.a+1)/(opts.par.b);
v_star = u_star^2;


J = gm_jac(u_star, v_star, opts.par.b);
opts.par.gamma = getCritLam(opts.par.d, pi*pi, J);

opts.par.a_heterogeneity = ...
    @(theta,x) (1 + sqrt(eps(1))*theta *cos(2*pi .* (x+0.5)));

opts.L = 0.5;
opts.Nx = 1e3;

p = init_gm_orig(opts);

% Set continuation parameters

p.nc.ilam=1;
p.nc.dsmin = 1e-5;
p.sol.ds=1e-6;
p.nc.dsmax=1e-1;
p.nc.dlammin = 1e-5;
p.nc.dlammax = 1e-1;
p.nc.lammin=0;
p.nc.lammax=1;

%p.sw.eigssol = 0;
p.sw.bifcheck = 1;
p.sw.spcalc = 1;
p.nc.tol = 1e-8;
p.fsol.tol = 1e-10;
p.fsol.fsol = 0;

p.fsol.disp = 0;

p = cont(p, 1000);



