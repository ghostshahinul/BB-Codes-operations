function edges = normalizeEdgesStruct(edges)
% normalizeEdgesStruct  Ensure cellulation bookkeeping fields exist.
    E = edges.numEdges;
    if ~isfield(edges, 'isCellulation') || isempty(edges.isCellulation)
        edges.isCellulation = false(E,1);
    end
    if ~isfield(edges, 'support')
        edges.support = [];
    end
end
