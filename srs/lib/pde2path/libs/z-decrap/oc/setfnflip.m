function fn=setfnflip(sd0,sp0,sd1,sp1,flip)
% setfnflip: convenience function to fill fn with dir and ptnames of CSS
% if flip=1 then flip start and target
%
%  fn=setfnflip(sd0,sp0,sd1,sp1,flip)
if flip==1; dud=sd0; dup=sp0; sd0=sd1; sp0=sp1; sd1=dud; sp1=dup; end % FLIP 
fn.sd0=sd0; fn.sp0=sp0; fn.sd1=sd1; fn.sp1=sp1; 