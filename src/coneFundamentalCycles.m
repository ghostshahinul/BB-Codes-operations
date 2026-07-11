function cycles = coneFundamentalCycles(D1)
% coneFundamentalCycles
% Fundamental cycle basis of the undirected graph represented by signed D1.
% Rows of D1 are edges; columns are vertices.

    D1 = double(D1);
    E = size(D1,1);
    V = size(D1,2);

    parentUF = 1:V;
    treeEdges = [];
    nonTreeEdges = [];

    for e = 1:E
        [u, v] = edgeEndpointsFromD1(D1, e);
        ru = ufFind(parentUF, u);
        rv = ufFind(parentUF, v);
        if ru ~= rv
            parentUF = ufUnion(parentUF, ru, rv);
            treeEdges(end+1) = e; %#ok<AGROW>
        else
            nonTreeEdges(end+1) = e; %#ok<AGROW>
        end
    end

    adj = cell(1,V);
    for idx = 1:numel(treeEdges)
        e = treeEdges(idx);
        [u, v] = edgeEndpointsFromD1(D1, e);
        adj{u} = [adj{u}; v, e];
        adj{v} = [adj{v}; u, e];
    end

    cycles = struct('vertices',{},'edges',{});

    for idx = 1:numel(nonTreeEdges)
        e = nonTreeEdges(idx);
        [u, v] = edgeEndpointsFromD1(D1, e);
        [pathVerts, pathEdges] = treePathBetweenVertices(adj, u, v);

        cycles(end+1).vertices = pathVerts; %#ok<AGROW>
        cycles(end).edges = [pathEdges, e];
    end
end
