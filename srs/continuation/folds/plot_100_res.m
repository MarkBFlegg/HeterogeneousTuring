h = @(k_sqr, d, jac)d * k_sqr.^2 -(d * jac(1,1) +jac(2,2)).*k_sqr + det(jac);
alpha = 1;
jac = @(u1,u2, alpha) [
    (-u2.^2) (-2*u1.*u2);
    (u2.^2) (2* u1.*u2 - alpha)
];

re_lam = @(k_sqr, d,jac) 0.5*(...
    -(k_sqr.*(1+d) - trace(jac)) + ...
    abs(real(sqrt( ...
        (k_sqr .*(1+d) - trace(jac)).^2 ...
        - 4 * h(k_sqr,d,jac)...
    )))...
);

min_ksqr = @(d,jac) ( ...
    (d*jac(1,1) + jac(2,2)) ...
    - real(sqrt((d*jac(1,1) + jac(2,2)).^2 - 4*d*det(jac)))...
) / (2 * d);

max_ksqr = @(d,jac) ( ...
    (d*jac(1,1) + jac(2,2)) ...
    + real(sqrt((d*jac(1,1) + jac(2,2)).^2 - 4*d*det(jac)))...
) / (2 * d);


load('100_res.mat');

min_n_vals = zeros(size(epsilon_vals));
max_n_vals = zeros(size(epsilon_vals));
gamma_vals = 1./(epsilon_vals.^2);

for i = 1:numel(epsilon_vals)
    min_n_vals(i) = sqrt(min_ksqr(1/40, gamma_vals(i) * jac(0.8, 1, 1)))/pi;
    max_n_vals(i) = sqrt(max_ksqr(1/40, gamma_vals(i) * jac(0.8, 1, 1)))/pi;
end

hold on
imagesc(epsilon_vals, n_vals, fold_lengths);
colormap gray
colorbar
plot(epsilon_vals, max_n_vals, 'Color','red','LineWidth',3);
plot(epsilon_vals, min_n_vals, 'Color','red','LineWidth',3);
xlim([min(epsilon_vals) max(epsilon_vals)])
ylim([min(n_vals) max(n_vals)])
title('Maximum $\theta$ before fold', 'Interpreter','latex', 'FontSize',20);
ylabel('n', 'Interpreter','latex', 'FontSize',20);
xlabel('$\epsilon$', 'Interpreter','latex', 'FontSize',20);
legend('Region where perturbation is unstable', 'Interpreter','latex', 'FontSize',14)
print('fold_size_n_epsilon.png', '-dpng', '-r400')