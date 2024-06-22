function [gamma] = solveCritGam(p)
%SOLVECRITEPS Summary of this function goes here
%   Detailed explanation goes here

    p = cont(p);
    fname = p.file.dir;
    pts = sort(getlabs(fname));
    p = loadp(fname, sprintf('pt%d', pts(end-1)));

    stable = min(real(p.sol.muv)) > 0;

    p.nc.ilam=2;
    p.nc.dsmin = 1e-10;
    p.nc.dsmax=1e-1;
    p.nc.dlammin = 1e-10;
    p.nc.dlammax = 1e-1;
    p.nc.lammin=0;
    p.nc.lammax=1000;
    
    p.myParams.stable = stable;
    if stable
        p.sol.ds=1e-10;
    else
        p.sol.ds=-1e-10;
    end
    p.fuha.ufu = @bifufu;
    close all
    p = cont(p);
    if p.file.bcount == 1
        %p = loadp(p.file.dir);
        gamma = nan; %getlam(p);
    else
        p = loadp(p.file.dir, sprintf('bpt%d', p.file.bcount-1));
        gamma = getlam(p);
    end
end

