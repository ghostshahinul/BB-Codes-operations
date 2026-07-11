function D1 = coneD1_signedOverZ(HZ, zi, edges)
% coneD1_signedOverZ  Signed incidence D1, with mapping-cone minus absorbed.
    HZ = double(HZ);
    [support, ~] = coneSupportZ(HZ, zi);
    w = numel(support);
    E = edges.numEdges;

    rows = [];
    cols = [];
    vals = [];

    for e = 1:E
        rows(end+1,1) = e;
        cols(end+1,1) = edges.uLocal(e);
        vals(end+1,1) = -edges.amount(e);

        rows(end+1,1) = e;
        cols(end+1,1) = edges.vLocal(e);
        vals(end+1,1) = +edges.amount(e);
    end

    D1 = sparse(rows, cols, vals, E, w);
end
