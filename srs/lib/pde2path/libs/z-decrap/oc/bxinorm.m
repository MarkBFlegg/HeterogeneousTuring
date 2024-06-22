function no=bxinorm(u,xi) 
% bxinorm: weighted norm (for iscarc)
n=length(u)-1; un=norm(u(1:n)); no=sqrt(xi*un^2+(1-xi)*u(n+1)^2);