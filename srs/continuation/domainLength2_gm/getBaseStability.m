function [stab] = getBaseStability(pp)
%GETBASESTABILITY Summary of this function goes here
%   Detailed explanation goes here
pp.file.dir
pp = cont(pp);

fname = pp.file.dir;
pts = sort(getlabs(fname));
pp = loadp(fname, sprintf('pt%d', pts(end-1)));
stab = min(real(pp.sol.muv));

end

