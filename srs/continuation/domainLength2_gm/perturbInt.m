function [y] = perturbInt(p)
%PERTURBINT Summary of this function goes here
%   Detailed explanation goes here
opts2 = odeset('Stats', 'on', 'OutputFcn', @dFunc, 'JPattern', sGjac_schk(p, p.u));
[~,y] = ode15s(@(t,y) -sG_schk(p, au2u(p, y)), [0 1e+12], u2au(p) + 0.01*rand(size(u2au(p))), opts2);
y = y(end, :);
end

