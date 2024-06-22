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

res = cell(3,1);

gamma = 3;
ksqr_vals = 0:0.01:40;
lam_vals = zeros(size(ksqr_vals));
for i = 1:numel(ksqr_vals)
    lam_vals(i) = re_lam(ksqr_vals(i), 1/40, gamma * jac(0.69, 1,1));
end
res{1}.gamma = gamma;
res{1}.ksqr_vals = ksqr_vals;
res{1}.lam_vals = lam_vals;

hold on
yline(0, 'LineWidth',2, 'Color', 'black');
plot(ksqr_vals, lam_vals, 'LineWidth',3)
xline(pi*pi, 'Color', 'red', 'LineWidth',3);
xline(pi*pi*2.^2, 'Color', 'red', 'LineWidth',3);
xlim([0 40])
ylim([-0.15 0.15])
xlabel("$k^2$", 'Interpreter','latex', 'FontSize',20);
ylabel("$\max_{i} \Re(\lambda_{k,i})$", 'Interpreter','latex', 'FontSize',20);
title("Dispersion relation with $\gamma = 3$", 'Interpreter','latex', 'FontSize',20);
legend('','','Eigenvalues', 'Interpreter','latex', 'FontSize',20)

figure(2)
gamma = 0.5;
ksqr_vals = 0:0.01:40;
lam_vals = zeros(size(ksqr_vals));
for i = 1:numel(ksqr_vals)
    lam_vals(i) = re_lam(ksqr_vals(i), 1/40, gamma * jac(0.69, 1,1));
end
res{2}.gamma = gamma;
res{2}.ksqr_vals = ksqr_vals;
res{2}.lam_vals = lam_vals;

hold on
yline(0, 'LineWidth',2, 'Color', 'black');
plot(ksqr_vals, lam_vals, 'LineWidth',3)
xline(pi*pi, 'Color', 'red', 'LineWidth',3);
xline(pi*pi*2.^2, 'Color', 'red', 'LineWidth',3);
xlim([0 40])
ylim([-0.15 0.15])
xlabel("$k^2$", 'Interpreter','latex', 'FontSize',20);
ylabel("$\max_{i} \Re(\lambda_{k,i})$", 'Interpreter','latex', 'FontSize',20);
title("Dispersion relation with $\gamma = 0.5$", 'Interpreter','latex', 'FontSize',20);
legend('','','Eigenvalues', 'Interpreter','latex', 'FontSize',20)



figure(3)
gamma = 1000;
ksqr_vals = 0:1:40000;
lam_vals = zeros(size(ksqr_vals));
for i = 1:numel(ksqr_vals)
    lam_vals(i) = re_lam(ksqr_vals(i), 1/40, gamma * jac(0.69, 1,1));
end
res{3}.gamma = gamma;
res{3}.ksqr_vals = ksqr_vals;
res{3}.lam_vals = lam_vals;

hold on
yline(0, 'LineWidth',2, 'Color', 'black');
plot(ksqr_vals, lam_vals, 'LineWidth',3)
n = 1;
while pi*pi*n*n <= max(ksqr_vals)
    xline(pi*pi*n*n, 'Color', 'red', 'LineWidth',1);
    n = n + 1;
end

xlim([0 2e4]);
ylim([-100 100])
xlabel("$k^2$", 'Interpreter','latex', 'FontSize',20);
ylabel("$\max_{i} \Re(\lambda_{k,i})$", 'Interpreter','latex', 'FontSize',20);
title("Dispersion relation with $\gamma = 1000$", 'Interpreter','latex', 'FontSize',20);
legend('','','Eigenvalues', 'Interpreter','latex', 'FontSize',20)

json_string = jsonencode(res);
disp_file = fopen('disp_plotting.json', 'w');
fprintf(disp_file, json_string);
fclose(disp_file);

