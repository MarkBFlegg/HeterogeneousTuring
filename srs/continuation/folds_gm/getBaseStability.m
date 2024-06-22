function [stab] = getBaseStability(pp)
%GETBASESTABILITY Summary of this function goes here
%   Detailed explanation goes here
pp.file.dir
pp = cont(pp);

fname = pp.file.dir;
pts = sort(getlabs(fname));
pp = loadp(fname, sprintf('pt%d', pts(end-1)));
final_theta = getlam(pp);
if final_theta < 0.99
    stab = "FOLD";
    return;
end
stab = min(real(pp.sol.muv));

end

