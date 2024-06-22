function [x,p] = lssGmres(A,b,p)
%LSSBICGSTAB Summary of this function goes here
%   Detailed explanation goes here
x = gmres(A,b);
end

