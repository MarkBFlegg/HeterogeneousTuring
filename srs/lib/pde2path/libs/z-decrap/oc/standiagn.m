function zdia=standiagn(sol,wnr) 
% standiagn: plot u_1, lam_1 for OC; use as template to plot adapted diagnostics 
global u0 u1 s1; 
x=sol.x; y=sol.y; sl=length(x); np=s1.np; zdia=zeros(2,sl); 
for i=1:sl
    zdia(1,i)=norm(y(1:np,i)-u0(1:np),'inf'); % diff state from IC 
    zdia(2,i)=norm(y(1:np,i)-u1(1:np),'inf'); % diff state from target
    zdia(3,i)=norm(y(np+1:2*np,i)-u0(np+1:2*np),'inf'); % diff costate from IC
    zdia(4,i)=norm(y(np+1:2*np,i)-u1(np+1:2*np),'inf'); % diff costate from target
end
figure(wnr); set(gca,'FontSize',s1.plot.fs); 
plot(x,zdia(1,:),'k',x,zdia(2,:),'b',x,zdia(3,:),'m',x,zdia(4,:),'g','LineWidth',2); 
legend('||u-u_0||_\infty','||u-u_1||_\infty','||\lambda-\lambda_0||_\infty','||\lambda-\lambda_1||_\infty'); 
axis tight;