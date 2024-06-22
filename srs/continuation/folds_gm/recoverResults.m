folders = dir('secant/crit_eps_*');

b_vals = linspace(0, 1 - 1.4901e-08, 1000);

results = {};
eps_vals = {};

for i = 1:numel(folders)
    names = split(folders(i).name, '_')
    if numel(names) <= 3
        continue;
    end
    
    results{str2num(names{3})}{str2num(names{4})+1} = folders(i);

    p = loadp(['secant/' folders(i).name]);
    eps_vals{str2num(names{3})}{str2num(names{4})+1} = p.u(p.nu+2);
end

crit_eps_vals_recovered = [];
crit_eps_vals_reference = [];
for i = 1:numel(eps_vals)
    if numel(eps_vals{i}) == 0
        crit_eps_vals_recovered(i) = -1;
        crit_eps_vals_reference(i) = -1;
        continue;
    end
    crit_eps_vals_recovered(i) = eps_vals{i}{end};
    crit_eps_vals_reference(i) = getCritEps(1, b_vals(i), 1/40, pi*pi);

end


