e = 1;
beta0 = 0.1;
p = makeP( ...
    sprintf('foldLength_analytical'), ...
    e, beta0...
);

function [p]=makeP(fname, epsilon, beta0)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.alpha = 1;
    opts.par.d = 1/40;
    opts.par.epsilon = epsilon;
    lScaling = 1;
    opts.par.betaf = @(x,theta) beta0 + (1- theta *cos(1*pi .* (x./lScaling+0.5)))/25;
    opts.par.eta = @(x,theta) 1 - opts.par.betaf(x,theta);
    opts.L = 0.5 *lScaling;
    opts.Nx = 1e4;
    
    p = init_schk(opts);
    p.myParams.lScaling = lScaling;
    
    % Set continuation parameters
    
    p.nc.ilam=1;
    p.nc.dsmin = 1e-5;
    p.sol.ds=1e-3;
    p.nc.dsmax=1e-1;
    p.nc.dlammin = 1e-5;
    p.nc.dlammax = 1e-1;
    p.nc.lammin=0;
    p.nc.lammax=1;
    
    p.sw.eigssol = 0;
    p.sw.bifcheck = 0;
    p.sw.spcalc = 1;
    p.nc.tol = 1e-8;
    p.fsol.tol = 1e-10;
    p.fsol.fsol = 0;
    
    p.fsol.disp = 0;
end