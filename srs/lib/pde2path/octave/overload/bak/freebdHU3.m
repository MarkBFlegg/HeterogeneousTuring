function es=freebdHU3(x,t) % emulates delaunayTriangulation.freeBoundary;  very inefficient! 
nt=size(t,1); e=[]; 
for i=1:nt
  e1=t(i,1:2); e2=t(i,2:3); e3=[t(i,3),t(i,1)]; % the three edges in row i
  bd1=1; bd2=1; bd3=1; % now check if these edges occur again ...
 for j=1:nt
   if j~=i; 
   f1=t(j,1:2); f2=t(j,2:3); f3=[t(j,3),t(j,1)];
   if (eqs(e1,f1) | eqs(e1,f2) | eqs(e1,f3)); bd1=0; end
   if (eqs(e2,f1) | eqs(e2,f2) | eqs(e2,f3)); bd2=0; end
   if (eqs(e3,f1) | eqs(e3,f2) | eqs(e3,f3)); bd3=0; end
 end 
 end
 if bd1; e=[e; e1]; end 
 if bd2; e=[e; e2]; end 
 if bd3; e=[e; e3]; end 
end
es=e; 
return 
e1=e, ne=size(e,1), % now resort such that each node appears once on each side 
es=e(1,:); ec=2:ne; % comparison rows 
for i=2:ne
  no=es(i-1,2); % last node 
  [r,c]=find(no==e(ec,1:2));  % find node in remaining rows 
  if c==1; es=[es; e(ec(r),:)];
  else es=[es; [e(ec(r),2) e(ec(r),1)]]; 
  end
  ec=setdiff(ec,i-1+r); 
  es, ec, pause 
end