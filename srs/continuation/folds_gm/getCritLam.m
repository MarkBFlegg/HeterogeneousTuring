function [lam]=getCritLam(d, k_sqr, J)
lam = (2*d*k_sqr) / ...
    ( (d*J(1,1) + J(2,2)) + sqrt((d*J(1,1) + J(2,2))^2 - 4*d*det(J)));
end