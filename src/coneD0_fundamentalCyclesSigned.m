function D0 = coneD0_fundamentalCyclesSigned(D1)
% coneD0_fundamentalCyclesSigned
% Signed fundamental cycle matrix before cellulation.

    cycles = coneFundamentalCycles(D1);
    E = size(D1,1);
    D0dense = zeros(numel(cycles), E);

    for c = 1:numel(cycles)
        D0dense(c,:) = signedCycleRowFromOrderedCycle(D1, cycles(c).vertices, cycles(c).edges);
    end

    D0 = sparse(D0dense);
end
