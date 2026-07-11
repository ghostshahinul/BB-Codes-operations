function [pathVerts, pathEdges] = treePathBetweenVertices(adj, startV, targetV)
% treePathBetweenVertices  BFS path in a tree/forest adjacency list.

    V = numel(adj);
    visited = false(1,V);
    prevV = zeros(1,V);
    prevE = zeros(1,V);

    q = startV;
    visited(startV) = true;

    while ~isempty(q)
        u = q(1);
        q(1) = [];
        if u == targetV
            break;
        end
        nbrs = adj{u};
        for j = 1:size(nbrs,1)
            v = nbrs(j,1);
            e = nbrs(j,2);
            if ~visited(v)
                visited(v) = true;
                prevV(v) = u;
                prevE(v) = e;
                q(end+1) = v; %#ok<AGROW>
            end
        end
    end

    if ~visited(targetV)
        error('No tree path found between vertices %d and %d.', startV, targetV);
    end

    vertsRev = targetV;
    edgesRev = [];
    cur = targetV;
    while cur ~= startV
        edgesRev(end+1) = prevE(cur); %#ok<AGROW>
        cur = prevV(cur);
        vertsRev(end+1) = cur; %#ok<AGROW>
    end

    pathVerts = fliplr(vertsRev);
    pathEdges = fliplr(edgesRev);
end
