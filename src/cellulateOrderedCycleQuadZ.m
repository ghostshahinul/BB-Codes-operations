function [D1_now, D0_now, F0_now, edges_now, addedList] = cellulateOrderedCycleQuadZ(D1_now, D0_now, F0_now, edges_now, vertices, edgeIDs, maxCycleLen)
% cellulateOrderedCycleQuadZ
% Split one ordered cycle into cycles of length <= maxCycleLen.
%
% For maxCycleLen=4, this repeatedly adds a chord from vertex 4 back to vertex 1:
%   [v1 v2 v3 v4 v5 ... vL]
% becomes
%   [v1 v2 v3 v4] plus [v1 v4 v5 ... vL].
%
% For a length-6 cycle this adds one chord and gives two 4-cycles, exactly the
% style needed for Example 11-type behavior.

    addedList = struct('edgeIndex',{},'tailLocal',{},'headLocal',{},'tailGlobal',{},'headGlobal',{});

    currVerts = vertices(:).';
    currEdges = edgeIDs(:).';

    while numel(currEdges) > maxCycleLen
        splitIdx = maxCycleLen;

        v1 = currVerts(1);
        vk = currVerts(splitIdx);

        % Add chord vk -> v1. This closes the first small cycle in the forward direction.
        [D1_now, F0_now, edges_now, newEdgeIdx] = addCellulationEdgeZ(D1_now, F0_now, edges_now, vk, v1);
        D0_now(:, end+1) = 0;  % pad previous D0 rows for the new edge

        % First small cycle: v1 -> v2 -> ... -> vk -> v1.
        smallVerts = currVerts(1:splitIdx);
        smallEdges = [currEdges(1:splitIdx-1), newEdgeIdx];
        row = signedCycleRowFromOrderedCycle(D1_now, smallVerts, smallEdges);
        D0_now(end+1, :) = row;

        % Remaining cycle: v1 -> vk -> v_{k+1} -> ... -> vL -> v1.
        % The first edge is the chord traversed opposite to its orientation.
        currVerts = [currVerts(1), currVerts(splitIdx:end)];
        currEdges = [newEdgeIdx, currEdges(splitIdx:end)];

        added.edgeIndex = newEdgeIdx;
        added.tailLocal = vk;
        added.headLocal = v1;
        if isfield(edges_now, 'support') && ~isempty(edges_now.support)
            added.tailGlobal = edges_now.support(vk);
            added.headGlobal = edges_now.support(v1);
        else
            added.tailGlobal = NaN;
            added.headGlobal = NaN;
        end
        addedList(end+1) = added; %#ok<AGROW>
    end

    if ~isempty(currEdges)
        row = signedCycleRowFromOrderedCycle(D1_now, currVerts, currEdges);
        D0_now(end+1, :) = row;
    end
end
