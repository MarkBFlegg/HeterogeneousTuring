close all; clear all;

b_vals = linspace(0, 24/25, 100);

crit_eps_vals = zeros(size(b_vals));

for betaI = 1:numel(b_vals)
    beta0 = b_vals(betaI);
    p = makeP(beta0, 1e-3, sprintf('crit_eps_%i', betaI));
    
    stab = getBaseStability(p)
    if stab > 0
        crit_eps_vals(betaI) = -1;
        continue;
    end
    p = makeP(beta0, 1, sprintf('crit_eps_%i', betaI));
    crit_eps_vals(betaI) = solveCritEps(p, 1);
end

function p = makeP(beta0, eps, fname)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.alpha = 1;
    opts.par.d = 1/40;
    opts.par.epsilon = eps;
    opts.par.betaf = @(x,theta) beta0 + (1 - theta *cos(pi .* (x+0.5)))/25;
    opts.par.eta = @(x,theta) 1 - opts.par.betaf(x,theta);
    opts.L = 0.5;
    opts.Nx = 1e4;
    
    p = init_schk(opts);
    p.myParams.lScaling = 1;
    
    % Set continuation parameters
    
    p.nc.ilam=1;
    p.nc.dsmin = 1e-5;
    p.sol.ds=1e-3;
    p.nc.dsmax=1e-1;
    p.nc.dlammin = 1e-5;
    p.nc.dlammax = 1e-1;
    p.nc.lammin=0;
    p.nc.lammax=1;
    
    p.usrlam = [1];
    
    p.sw.eigssol = 0;
    p.sw.bifcheck = 1;
    p.sw.spcalc = 1;
    p.nc.tol = 1e-8;
    p.fsol.tol = 1e-10;
    p.fsol.fsol = 0;
    
    p.fsol.disp = 0;

end
