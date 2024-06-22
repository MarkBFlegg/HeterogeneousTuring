function [x] = SolvePoorlyConditionedSystemDirect(A, b)
%SOLVEPOORLYCONDITIONEDSYSTEM Summary of this function goes here
%   Detailed explanation goes here

[P,R,C] = equilibrate(A);
B = R*P*A*C;
d = R*P*b;
y = B \ d;
x = C * y;

end

