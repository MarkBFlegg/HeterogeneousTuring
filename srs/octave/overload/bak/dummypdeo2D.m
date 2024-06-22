classdef dummypdeo2D < pde  
methods(Access = public)
  function o=stanpdeo2D(lx,ly,varargin) % constructor 
     o.grid=grid2D; o.grid.mySquare(lx,ly,varargin{1}); 
     o.fem=lagrange12D; 111
   end
end
methods(Access = protected) % only here since required by pde-class
    function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
end
end
