 % 
 x=0:0.25:1; y=0:0.5:1; [X,Y]=meshgrid(x,y); Z=0*X; 
  [V F] = surfToMesh(X, Y); 
 t=triangulateFaces(F); 
 figure(1); clf;   drawMesh(V, t);
 ec = meshEdges(V, t);
 inds = meshBoundaryEdgeIndices(V, ec, t);
 e=ec(inds,:); 
 ep = V(e(:,1), :);  
 figure(2); clf; plot(ep', '*', 'color', 'b');