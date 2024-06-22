close all; clear all;

b_vals = linspace(0.67, 1.0, 200);

crit_eps_vals = zeros(size(b_vals));
ref_eps_vals = zeros(size(b_vals));
ref_eps_vals_2 = zeros(size(b_vals));


for bI = 1:numel(b_vals)-1
    b0 = b_vals(bI);

    epsilon = getCritEps(1.0, b0, 1/40, pi*pi);
    if imag(epsilon) == 0
        ref_eps_vals(bI) = epsilon;
        ref_eps_vals_2(bI) = getCritEps(1.0, b0+sqrt(eps(1)), 1/40, pi*pi);
    end
    
    dirname = sprintf('crit_eps_%i', bI);
    if isfolder(dirname) && isfile([dirname '/bpt1.mat'])
        p = loadp(dirname, 'bpt1');
        crit_eps_vals(bI) = getlam(p);
        continue;
    end

    p = makeP(b0, 1e-2, dirname);
    
    stab = getBaseStability(p)
    if stab > 0
        crit_eps_vals(bI) = -1;
        continue;
    end


    p = makeP(b0, epsilon, sprintf('crit_eps_%i', bI));
    crit_eps_vals(bI) = solveCritEps(p);
end


function p = makeP(beta0, epsilon, fname)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.alpha = 1;
    opts.par.d = 1/40;
    opts.par.epsilon = epsilon;
    opts.par.betaf = @(x,theta) beta0 * (1 + theta*1.4901e-08 *cos(2*pi .* (x+0.5)));
    opts.par.eta = @(x,theta) 1 - opts.par.betaf(x,theta);
    opts.L = 0.5;
    opts.Nx = 1e2;
    
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
    
    %p.sw.eigssol = 0;
    p.sw.bifcheck = 1;
    p.sw.spcalc = 1;

    x = getpte(p).';
    p.nc.tol = 1e-12/(x(2)-x(1));
    p.fsol.tol = 1e-16;
    p.fsol.fsol = 0;
    
    p.fsol.disp = 0;

end
