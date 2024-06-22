function dt = dt2(X,Y)
   [V F] = surfToMesh(X, Y); 
 t=triangulateFaces(F); 
 ec = meshEdges(V, t);
 inds = meshBoundaryEdgeIndices(V, ec, t);
 e=ec(inds,:); 
 dt.ConnectivityList=t; dt.freeBoundary=e; %dt.Points=x; 
%t=delaunayn(x); t
%e=freebdHU(x,t); 
dt.ConnectivityList=t; dt.freeBoundary=e; %dt.Points=x; 
  end
