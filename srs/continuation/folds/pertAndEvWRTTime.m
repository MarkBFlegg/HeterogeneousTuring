function q = pertAndEvWRTTime(p)
    q = p;
    opts2 = odeset('Stats', 'on', 'OutputFcn', @dFunc, 'JPattern', sGjac_schk(p, p.u));
    [~,y] = ode15s(@(t,y) -sG_schk(p, au2u(p, y)), [0 1e+12], u2au(p) + 0.01*rand(size(u2au(p))), opts2);
    q.u = au2u(p,y(end, :));

end