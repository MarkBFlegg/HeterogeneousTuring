function [eps] = solveCritEps(p, start)
%SOLVECRITEPS Summary of this function goes here
%   Detailed explanation goes here
    
    counter = 0;
    eps = fsolve(@(x) solve_f(p, x), start);

    function stab=solve_f(p, test_eps) 
        q = p;
        q = setfn(q, sprintf('%s_%d', p.file.dir, counter));
        counter = counter + 1;
        q.u(q.nu+2) = test_eps;
        q.u(q.nu+1) = 0;
        stab = getBaseStability(q)
    end
end

