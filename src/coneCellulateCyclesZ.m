function [D1_cell, D0_cell, F0_cell, edges_cell, cellInfo] = coneCellulateCyclesZ(D1_pre, F0_pre, edges_pre, maxCycleLen)
% coneCellulateCyclesZ
% Find a fundamental cycle basis of the graph represented by D1_pre and
% split every cycle with length > maxCycleLen by adding chord edges.
%
% New chord edges:
%   - are appended as extra rows of D1
%   - have zero columns in F0
%   - are recorded in edges_cell as isCellulation=true
%
% D0_cell is rebuilt from the cellulated small cycles and satisfies D0_cell*D1_cell=0 over Z.

    if nargin < 4 || isempty(maxCycleLen)
        maxCycleLen = 4;
    end

    D1_cell = full(double(D1_pre));
    F0_cell = full(double(F0_pre));
    edges_cell = normalizeEdgesStruct(edges_pre);

    cycles = coneFundamentalCycles(D1_pre);
    D0_dense = zeros(0, size(D1_cell,1));

    cellInfo.originalCycles = cycles;
    cellInfo.addedEdges = [];
    cellInfo.numAddedEdges = 0;
    cellInfo.maxCycleLen = maxCycleLen;

    for c = 1:numel(cycles)
        verts = cycles(c).vertices;
        eids  = cycles(c).edges;

        [D1_cell, D0_dense, F0_cell, edges_cell, added] = ...
            cellulateOrderedCycleQuadZ(D1_cell, D0_dense, F0_cell, edges_cell, verts, eids, maxCycleLen);

        for a = 1:numel(added)
            cellInfo.numAddedEdges = cellInfo.numAddedEdges + 1;
            cellInfo.addedEdges = [cellInfo.addedEdges; added(a)]; %#ok<AGROW>
        end
    end

    D1_cell = sparse(D1_cell);
    D0_cell = sparse(D0_dense);
    F0_cell = sparse(F0_cell);

    coneCheckCyclesZ(D0_cell, D1_cell);
end
