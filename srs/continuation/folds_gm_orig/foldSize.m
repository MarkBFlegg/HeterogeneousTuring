close all; clear all;
length_vals = linspace(2,20,100);
a0_vals = linspace(0.0, 0.4, 50);

[ll, aa] = meshgrid(length_vals, a0_vals);
fold_lengths = zeros(size(ll));

parfor i = 1:numel(ll)
    fname = sprintf('foldLength1_%d', i);

    % uncomment for run 2.
    %if isfolder(fname) && isfile([fname '/fpt1.mat'])
    %    pp = loadp(fname, 'fpt1');
    %    if getlam(pp) < 1.0
    %        rmdir(fname, 's');
    %    end
    %end

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

                pp = cont(pp);
                pp.nc.dlammax = 1;
                pp.nc.dsmax = 1;
                l=getlam(pp);
            end
        end
        fold_lengths(i) = l;
        continue;
    end
    a0 = aa(i);
    l = ll(i);
    p = makeP( ...
        sprintf('foldLength1_%d', i), ...
        l, a0...
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

function [p]=makeP(fname, l, a0)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.a = a0;
    opts.par.b = 1;
    opts.par.d = 20;
    opts.par.gamma = l*l;
    opts.par.a_heterogeneity = @(theta,x) (1+ theta *cos(pi .* (x+0.5)));
    opts.L = 0.5;
    opts.Nx = 2e2;
    
    p = init_gm_orig(opts);
    
    % Set continuation parameters
    
    p.nc.ilam=1;
    p.nc.dsmin = 1e-10;
    p.sol.ds=1e-10;
    p.nc.dsmax=1e-1;
    p.nc.dlammin = 1e-10;
    p.nc.dlammax = 1e-1;
    p.nc.lammin=0;
    p.nc.lammax=20;

    % Run #2 was run with 
    % p.nc.dsmax =1e-2;
    % p.nc.dlammax = 1e-2;
    
    p.sw.eigssol = 0;
    p.sw.bifcheck = 0;
    p.sw.spcalc = 1;
    p.nc.tol = 1e-8;
    p.fsol.tol = 1e-10;
    p.fsol.fsol = 0;
    
    p.fsol.disp = 0;
    p.fuha.ufu = @foldufu;
end