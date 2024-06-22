function [eps] = getCritEps(alpha, beta0, d, k_sqr)
%GETCRITEPS Summary of this function goes here
%   Detailed explanation goes here
eps = sqrt(1/getCritLam(d, k_sqr, jac(alpha, beta0*alpha*alpha, 1/alpha)));


    function j=jac(alpha, u1,u2)
        j = [-u2*u2 -2*u1*u2; u2*u2 2*u1*u2-alpha];
    end

end

