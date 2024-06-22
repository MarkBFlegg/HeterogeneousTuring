close all; clear all;

b_vals = linspace(0.67, 1.0, 100);
theta_vals = linspace(-0.5, 0.5, 7);

[bb, tt] = meshgrid(b_vals, theta_vals);

crit_eps_vals = zeros(size(bb));
ref_eps_vals = zeros(size(bb));
ref_eps_vals_2 = zeros(size(bb));


for bI = 201:201 %1:numel(bb)
    b0 = bb(bI);
    theta = tt(bI);

    epsilon = getCritEps(1.0, b0, 1/40, pi*pi);
    ref_eps_vals(bI) = epsilon;
    
    dirname = sprintf('crit_eps_%ib', bI);
    if isfolder(dirname) && isfile([dirname '/bpt2.mat'])
        p = loadp(dirname, 'bpt2');
        crit_eps_vals(bI) = getlam(p);
        continue;
    elseif isfolder(dirname) && isfile([dirname '/bpt1.mat'])
        p = loadp(dirname, 'bpt1');
        crit_eps_vals(bI) = getlam(p);
        continue;
    end

    p = makeP(b0, 1e-2, sprintf('crit_eps_%ia', bI), theta);
    
    stab = getBaseStability(p)
    if stab > 0
        crit_eps_vals(bI) = -1;
        continue;
    end
    if imag(epsilon) == 0
        eStart = epsilon;
        p = makeP(b0, epsilon, sprintf('crit_eps_%ib', bI), theta);
    else
        eStart = abs(epsilon);
        p = makeP(b0, real(epsilon), sprintf('crit_eps_%ib', bI), theta);
    end
    crit_eps_vals(bI) = solveCritEps(p);
    %crit_eps_vals(bI) = fzero(@(e) getBaseStability(makeP(b0, e, sprintf('crit_eps_%ia_%f', bI, e), theta)), ...
    %    eStart);
end


function p = makeP(beta0, epsilon, fname, theta_max)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.alpha = 1;
    opts.par.d = 1/40;
    opts.par.epsilon = epsilon;
    opts.par.betaf = @(x,theta) beta0 * (1 + theta_max*theta *cos(2*pi .* (x+0.5)));
    opts.par.eta = @(x,theta) 1 - opts.par.betaf(x,theta);
    opts.L = 0.5;
    opts.Nx = 1e3;
    
    p = init_schk(opts);
    p.myParams.lScaling = 1;
    
    % Set continuation parameters
    
    p.nc.ilam=1;
    p.nc.dsmin = 1e-5;
    p.sol.ds=1e-5;
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
