close all; clear all;

a_vals = linspace(0, 0.5, 200);

crit_gam_vals = zeros(size(a_vals));
reference_gam_vals = zeros(size(a_vals));
reference_gam_vals_2 = zeros(size(a_vals));

parfor aI = 2:numel(a_vals)
    a0 = a_vals(aI);
    
    ref_gamma = getCritLam(20, pi*pi, gm_jac(1+a0, (1+a0)^2, 1.0));
    ref_gamma_2 = getCritLam(20, pi*pi, gm_jac(1+a0+sqrt(eps(1)), (1+a0+sqrt(eps(1)))^2, 1.0));
    
    if imag(ref_gamma) == 0
        reference_gam_vals(aI) = ref_gamma;
        reference_gam_vals_2(aI) = ref_gamma_2;
    end

    dirname = sprintf('crit_eps_%i', aI);
    if isfolder(dirname) && isfile([dirname '/bpt1.mat'])
        p = loadp(dirname, 'bpt1');
        crit_gam_vals(aI) = getlam(p);
    end
    p = makeP(a0, 1000, dirname);
    
    stab = getBaseStability(p)
    if stab > -1e-5
        crit_gam_vals(aI) = -1;
        continue;
    end
    l = sqrt(reference_gam_vals(aI));

    p = makeP(a0, l, sprintf('crit_eps_%i', aI));
    crit_gam_vals(aI) = solveCritGam(p);
end

function p = makeP(a0, l, fname)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.a = a0;
    opts.par.b = 1.0;
    opts.par.d = 20;
    opts.par.gamma = l*l;
    opts.par.a_heterogeneity = @(theta,x) 1 + sqrt(eps(1)) * theta *cos(2*pi .* (x+0.5));
    opts.L = 0.5;
    opts.Nx = 2e2;
    
    p = init_gm_orig(opts);
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
