 % create centered icosahedron
     [v, f] = createIcosahedron;
     v(:,3) = v(:,3) - mean(v(:,3));
     % convert to simili-sphere
     [v2, f2] = subdivideMesh(v, f, 2);
     v3 = normalizeVector3d(v2);
     % clip with plane
    % f2=f; v3=v; 
     plane = createPlane([0 0 0], [-1 -2 3]);
     [vc, fc] = clipMeshVertices(v3, f2, plane, 'shape', 'plane');
   % vc=v; fc=f; 
     figure(1); clf; drawMesh(vc, fc); axis equal; view(3);
     % draw boundary vertices
     ec = meshEdges(vc, fc);
     inds = meshBoundaryEdgeIndices(vc, ec, fc);
     edges = [vc(ec(inds, 1), :) vc(ec(inds, 2), :)];
     hold on; drawEdge3d(edges, 'linewidth', 2, 'color', 'b');