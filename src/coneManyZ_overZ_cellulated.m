function [HX_new, HZ_new, allData] = coneManyZ_overZ_cellulated(HX, HZ, zRows, maxCycleLen)
% coneManyZ_overZ_cellulated
% Simultaneous signed/integer coning of many Z rows, INCLUDING CELLULATION.
%
% zRows refer to rows of the ORIGINAL HZ matrix.

    if nargin < 4 || isempty(maxCycleLen)
        maxCycleLen = 4;
    end

    HX = double(HX);
    HZ = double(HZ);
    zRows = unique(zRows(:).');

    checkIntegerMatrix(HX, 'HX');
    checkIntegerMatrix(HZ, 'HZ');
    checkCSSCommutationZ(HX, HZ, 'input HX,HZ');

    nX = size(HX,1);
    n  = size(HX,2);
    nZ = size(HZ,1);

    if any(zRows < 1) || any(zRows > nZ)
        error('zRows contains an invalid Z-row index.');
    end

    blocksD0 = {};
    F0_list = {};
    D1_list = {};
    F1_list = {};
    supportSizes = zeros(1,numel(zRows));
    edgeSizes = zeros(1,numel(zRows));
    cycleSizes = zeros(1,numel(zRows));
    allData = struct();

    for t = 1:numel(zRows)
        zi = zRows(t);
        [support, zSigns] = coneSupportZ(HZ, zi);
        F1 = coneF1_signed(HZ, zi);

        edges0 = coneEdges_signedOverZ(HX, HZ, zi);
        D1_pre = coneD1_signedOverZ(HZ, zi, edges0);
        F0_pre = coneF0_fromEdges(nX, edges0);
        D0_pre = coneD0_fundamentalCyclesSigned(D1_pre);

        coneCheckCoreZ(HX, F1, F0_pre, D1_pre);
        coneCheckCyclesZ(D0_pre, D1_pre);

        [D1, D0, F0, edges, cellInfo] = coneCellulateCyclesZ(D1_pre, F0_pre, edges0, maxCycleLen);

        coneCheckCoreZ(HX, F1, F0, D1);
        coneCheckCyclesZ(D0, D1);

        blocksD0{t} = D0;
        F0_list{t} = F0;
        D1_list{t} = D1;
        F1_list{t} = F1;
        supportSizes(t) = numel(support);
        edgeSizes(t) = size(D1,1);
        cycleSizes(t) = size(D0,1);

        allData(t).zi = zi;
        allData(t).support = support;
        allData(t).zSigns = zSigns;
        allData(t).F1 = F1;
        allData(t).F0 = F0;
        allData(t).D1 = D1;
        allData(t).D0 = D0;
        allData(t).D1_pre = D1_pre;
        allData(t).D0_pre = D0_pre;
        allData(t).F0_pre = F0_pre;
        allData(t).edges = edges;
        allData(t).edges_pre = edges0;
        allData(t).cellInfo = cellInfo;
        allData(t).numVertices_Q = numel(support);
        allData(t).numEdges_X = size(D1,1);
        allData(t).numCycles_R = size(D0,1);
        allData(t).numCellulationEdges = cellInfo.numAddedEdges;
    end

    Etotal = sum(edgeSizes);
    Ctotal = sum(cycleSizes);
    Wtotal = sum(supportSizes);

    D0_block = sparse(Ctotal, Etotal);
    F0_all = sparse(nX, Etotal);
    HZ_top = sparse(Wtotal, Etotal + n);

    edgeOffset = 0;
    cycleOffset = 0;
    zrowOffset = 0;

    for t = 1:numel(zRows)
        D0 = blocksD0{t};
        F0 = F0_list{t};
        D1 = D1_list{t};
        F1 = F1_list{t};

        E = size(D1,1);
        C = size(D0,1);
        W = size(D1,2);

        if C > 0
            D0_block(cycleOffset+(1:C), edgeOffset+(1:E)) = D0;
        end
        F0_all(:, edgeOffset+(1:E)) = F0;

        HZ_top(zrowOffset+(1:W), edgeOffset+(1:E)) = D1.';
        HZ_top(zrowOffset+(1:W), Etotal+(1:n)) = F1.';

        edgeOffset = edgeOffset + E;
        cycleOffset = cycleOffset + C;
        zrowOffset = zrowOffset + W;
    end

    HZ_remaining = HZ;
    HZ_remaining(zRows,:) = [];

    HX_new = [D0_block, sparse(Ctotal, n); F0_all, sparse(HX)];
    HZ_new = [HZ_top; sparse(size(HZ_remaining,1), Etotal), sparse(HZ_remaining)];

    checkCSSCommutationZ(HX_new, HZ_new, 'simultaneously coned HX_new,HZ_new');
end
