close all; clear all;
epsilon_vals = linspace(0.01, 2, 100);
beta0_vals = linspace(0.6, 1-1/25, 50);
n_vals = 1:10;
%fold_lengths = zeros(numel(beta0_vals), numel(epsilon_vals));

[ee, nn] = meshgrid(epsilon_vals, n_vals);
fold_lengths = zeros(size(ee));

for i = 101:101 %1:numel(ee)
    fname = sprintf('foldLength_%d', i)
    beta0 = 0.8;
    e = ee(i);
    n = nn(i);
    p = makeP( ...
        fname, ...
        e, beta0,n...
    );
    fold_lengths(i) = getFoldLength(p);
    %figure(3)
    %imagesc(fold_lengths);
end

function l=getFoldLength(pp)
    pp = cont(pp);

    fname = pp.file.dir;
    pts = sort(getlabs(fname));

    if numel(pts) < 2
        l = 0.0;
        return;
    end
    pp = loadp(fname, sprintf('pt%d', pts(end-1)));
    final_theta = getlam(pp);
    if pp.file.fcount > 1
        
        pp = loadp(fname, sprintf('fpt1'));
        l = getlam(pp);
        return;
    end
    l=final_theta;
end

function [p]=makeP(fname, epsilon, beta0,n)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.alpha = 1;
    opts.par.d = 1/40;
    opts.par.epsilon = epsilon;
    lScaling = 1;
    opts.par.betaf = @(x,theta) beta0 * (1+ theta *cos(n*pi .* (x./lScaling+0.5)));
    opts.par.eta = @(x,theta) 1 - opts.par.betaf(x,theta);
    opts.L = 0.5 *lScaling;
    opts.Nx = 1e3;
    
    p = init_schk(opts);
    p.myParams.lScaling = lScaling;
    
    % Set continuation parameters
    
    p.nc.ilam=1;
    p.nc.dsmin = 1e-5;
    p.sol.ds=1e-3;
    p.nc.dsmax=1e-2;
    p.nc.dlammin = 1e-5;
    p.nc.dlammax = 1e-2;
    p.nc.lammin=0;
    p.nc.lammax=1;
    
    p.sw.eigssol = 0;
    p.sw.bifcheck = 0;
    p.sw.spcalc = 1;
    p.nc.tol = 1e-8;
    p.fsol.tol = 1e-10;
    p.fsol.fsol = 0;
    
    p.fsol.disp = 0;
    p.fuha.ufu = @foldufu;
end