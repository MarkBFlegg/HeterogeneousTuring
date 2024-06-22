function [J] = gm_jac(u1, u2, b)
%GM_JAC Summary of this function goes here
%   Detailed explanation goes here

f1u = ( ...
    - b ...
    + 2 * u1 ./ u2 ...
);
f1v = ( ...
    - u1.^2 ./ (u2.^2) ...
);
f2u = ( ...
    2 * u1 ...
);
f2v = ( ...
    - 1 ...
);

J = [f1u f1v; f2u f2v];

end
