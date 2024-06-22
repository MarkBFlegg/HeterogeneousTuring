function [x,p] = lssBicgstab(A,b,p)
%LSSBICGSTAB Summary of this function goes here
%   Detailed explanation goes here
x = bicgstab(A,b);
end

