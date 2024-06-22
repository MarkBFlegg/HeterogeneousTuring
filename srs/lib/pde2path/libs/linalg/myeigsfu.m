function y=myeigsfu(p,A,B,sig,b) 
% gcafun: solver to be called in eigs; calls ilueigs 
A=A-sig*B; [y,~]=myIterativeLss(A,b,p); 