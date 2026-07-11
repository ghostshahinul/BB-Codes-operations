function F0 = coneF0_fromEdges(nX, edges)
% coneF0_fromEdges  Ownership map f0. ownerXRow=0 means zero column.
    E = edges.numEdges;
    if E == 0
        F0 = sparse(nX, 0);
        return;
    end
    valid = edges.ownerXRow > 0;
    F0 = sparse(edges.ownerXRow(valid), find(valid), 1, nX, E);
end
