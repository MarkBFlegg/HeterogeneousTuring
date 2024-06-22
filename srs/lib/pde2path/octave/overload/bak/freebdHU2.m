function e=freebdHU2(x,t)
  t, x
  if 1
alle=[t(:,[3 1]); t(:,[2,3]); t(:,[1 2])]; alle 
e=unique(alle,'rows'); e, pause 
return 
end 
nt=size(t,1); e=[]; 
for i=1:nt
  e1=t(i,1:2); e2=t(i,2:3); e3=[t(i,3),t(i,1)]; % the three edges in roow i
  bd1=1; bd2=1; bd3=1; % now check if these edges occur again ...
 for j=1:nt
   if j==i; j=j+1; end % don't compare with same row 
   try; %if j<nt+1
   f1=t(j,1:2); f2=t(j,2:3); f3=[t(j,3),t(j,1)];
   if (isequal(e1,f1) | isequal(e1,f2) | isequal(e1,f3)); bd1=0; end
   if (isequal(e2,f1) | isequal(e2,f2) | isequal(e2,f3)); bd2=0; end
   if (isequal(e3,f1) | isequal(e3,f2) | isequal(e3,f3)); bd3=0; end
   end
 end
 if bd3; e=[e; e3]; end 
 if bd2; e=[e; e2]; end 
 if bd1; e=[e; e1]; end 
end
end