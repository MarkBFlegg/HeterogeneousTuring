close all; clear all;
length_vals = linspace(1,200,100);
a0_vals = linspace(0.0, 1-1/25, 50);
n_vals = 1:10;
%fold_lengths = zeros(numel(beta0_vals), numel(epsilon_vals));

[ll, nn] = meshgrid(length_vals, n_vals);
fold_lengths = zeros(size(ll));

for i = 1:numel(ll)
    fname = sprintf('foldLength_%d', i);
    if isfolder(fname)
        pts = sort(getlabs(fname));
    
        if numel(pts) < 2
            l = 0.0;
        else
            pp = loadp(fname, sprintf('pt%d', pts(end-1)));
            final_theta = getlam(pp);
            if pp.file.fcount > 1
                pp = loadp(fname, sprintf('fpt1'));
                l = getlam(pp);
            else
                l=final_theta;
            end
        end
        fold_lengths(i) = l;
        continue;
    end
    a0 = 0.1;
    l = ll(i);
    n = nn(i);
    p = makeP( ...
        fname, ...
        l, a0,n...
    );
    fold_lengths(i) = getFoldLength(p);
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

function [p]=makeP(fname, l, a0,n)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.a = a0;
    opts.par.b = 1.0;
    opts.par.d = 20;
    opts.par.gamma = l*l;
    opts.par.a_heterogeneity = @(theta,x) (1+ theta *cos(n*pi .* (x+0.5)));
    opts.L = 0.5 ;
    opts.Nx = 1e3;
    
    p = init_gm_orig(opts);
    
    % Set continuation parameters
    
    p.nc.ilam=1;
    p.nc.dsmin = 1e-5;
    p.sol.ds=1e-10;
    p.nc.dsmax=1e-1;
    p.nc.dlammin = 1e-5;
    p.nc.dlammax = 1e-1;
    p.nc.lammin=0;
    p.nc.lammax=10;
    
    p.sw.eigssol = 0;
    p.sw.bifcheck = 0;
    p.sw.spcalc = 1;
    p.nc.tol = 1e-8;
    p.fsol.tol = 1e-10;
    p.fsol.fsol = 0;
    
    p.fsol.disp = 0;
    p.fuha.ufu = @foldufu;
end