% QUALITATIVE RESULTS
close all; clear all;

res = {};

%
L = 1;
epsilon = 1/L;
beta0 = 0.8;
n = 1;
p = makeP('qualitative_1', epsilon, beta0, n);
p = cont(p);

%p = loadp('qualitative_1', 'pt49');
%p.sol.ds = -p.sol.ds;
%p.u(1:p.np) = flipud(p.u(1:p.np));
%p.u(p.np+1:2*p.np) = flipud(p.u(p.np+1:2*p.np));
%p.nc.lammax = 0.1;
%p = cont(p);
%}


%{
L = 1.9;
epsilon = 1/L;
beta0 = 0.8;
n = 1;
p = makeP('qualitative_2', epsilon, beta0, n);

p.nc.lammin=-1;
p = cont(p, 280);
%}

%{
L = 3;
epsilon = 1/L;
beta0 = 0.8;
n = 1;
p = makeP('qualitative_3', epsilon, beta0, n);
p = cont(p);
%}


%{
L = 30;
epsilon = 1/L;
beta0 = 0.8;
n = 1;
p = makeP('qualitative_4', epsilon, beta0, n);
p = cont(p);
%}


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

function [p]=makeP(fname, epsilon, beta0, n)
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
    p.sol.ds=1e-2;
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
    p.file.smod = 10;
end