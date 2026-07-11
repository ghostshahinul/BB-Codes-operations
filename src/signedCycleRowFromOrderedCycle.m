function row = signedCycleRowFromOrderedCycle(D1, vertices, edgeIDs)
% signedCycleRowFromOrderedCycle
% Build one signed D0 row for an ordered cycle.
%
% vertices = [v1 ... vL]
% edgeIDs(k) connects vertices(k) to vertices(k+1), with vertices(L+1)=vertices(1).
%
% Coefficient rule:
%   D1 edge orientation tail->head is encoded as -1 at tail, +1 at head.
%   If traversal direction matches tail->head, coefficient is +1.
%   If traversal direction is opposite, coefficient is -1.

    D1 = double(D1);
    Etotal = size(D1,1);
    L = numel(edgeIDs);

    if numel(vertices) ~= L
        error('Number of vertices and edgeIDs must match for a closed cycle representation.');
    end

    row = zeros(1,Etotal);

    for k = 1:L
        a = vertices(k);
        if k < L
            b = vertices(k+1);
        else
            b = vertices(1);
        end

        e = edgeIDs(k);
        [tail, head, amount] = edgeOrientationFromD1(D1, e);
        if amount ~= 1
            error('Cellulation currently expects unit-weight signed edges. Edge %d has amount %g.', e, amount);
        end

        if tail == a && head == b
            coeff = +1;
        elseif tail == b && head == a
            coeff = -1;
        else
            error('Edge %d does not connect consecutive vertices %d and %d.', e, a, b);
        end

        row(e) = coeff;
    end

    test = row * D1;
    if any(test ~= 0)
        error('Constructed cycle row does not satisfy row*D1=0. max=%g', max(abs(test)));
    end
end
