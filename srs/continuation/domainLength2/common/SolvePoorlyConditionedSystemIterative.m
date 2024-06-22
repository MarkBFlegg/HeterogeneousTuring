function [x] = SolvePoorlyConditionedSystemIterative(A, b, tol, maxiter)
%SOLVEPOORLYCONDITIONEDSYSTEM Summary of this function goes here
%   Detailed explanation goes here

[P,R,C] = equilibrate(A);
B = R*P*A*C;
d = R*P*b;


%setup = struct('type','ilutp','droptol',1e-6);
%[L,U] = ilu(B,setup);
%[y,~,~,~,~] = bicgstabl(B,d,tol,maxiter,L,U);
[y,~,~,~,~] = bicgstabl(B,d,tol,maxiter);
x = C * y;


end

