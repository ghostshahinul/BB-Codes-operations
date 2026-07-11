function [HX_new, HZ_new, data] = coneOneZ_overZ_cellulated(HX, HZ, zi, maxCycleLen)
% coneOneZ_overZ_cellulated
% Signed/integer coning of one Z-check row over Z, INCLUDING CELLULATION.
%
% Inputs:
%   HX          : nX x n integer X-check matrix
%   HZ          : nZ x n integer/signed Z-check matrix, e.g. HZ_ldi in {-1,0,+1}
%   zi          : MATLAB 1-based row index of the Z-check to cone
%   maxCycleLen : optional maximum cycle/check weight for D0 rows; default 4
%
% Outputs:
%   HX_new : coned X-check matrix over Z
%   HZ_new : coned Z-check matrix over Z
%   data   : struct containing support, F1, F0, D1, D0, edges, cellulation info
%
% Sign convention:
%   The mapping-cone minus is already absorbed into D1.
%   Therefore HZ_new is assembled as [D1.' F1.']; do NOT use [-D1.' F1.'].
%
% Required exact checks over Z:
%   F0*D1 + HX*F1 == 0
%   D0*D1 == 0
%   HX_new*HZ_new.' == 0
%
% Cellulation convention:
%   Every cycle longer than maxCycleLen is split by adding chord edges.
%   New chord edges are added to Xi, hence new qubits, but they have no owner X-check,
%   so their columns in F0 are zero. This follows the paper's (_,j,k) convention.

    if nargin < 4 || isempty(maxCycleLen)
        maxCycleLen = 4;
    end
    if maxCycleLen < 3
        error('maxCycleLen must be at least 3. Use 4 for the paper-style cellulation target.');
    end

    HX = double(HX);
    HZ = double(HZ);

    checkIntegerMatrix(HX, 'HX');
    checkIntegerMatrix(HZ, 'HZ');
    checkCSSCommutationZ(HX, HZ, 'input HX,HZ');

    [support, zSigns] = coneSupportZ(HZ, zi);
    F1 = coneF1_signed(HZ, zi);

    % Pre-cellulation objects.
    edges0 = coneEdges_signedOverZ(HX, HZ, zi);
    D1_pre = coneD1_signedOverZ(HZ, zi, edges0);
    F0_pre = coneF0_fromEdges(size(HX,1), edges0);
    D0_pre = coneD0_fundamentalCyclesSigned(D1_pre);

    coneCheckCoreZ(HX, F1, F0_pre, D1_pre);
    coneCheckCyclesZ(D0_pre, D1_pre);

    % Cellulate long cycles and rebuild D0 for the cellulated graph.
    [D1, D0, F0, edges, cellInfo] = coneCellulateCyclesZ(D1_pre, F0_pre, edges0, maxCycleLen);

    coneCheckCoreZ(HX, F1, F0, D1);
    coneCheckCyclesZ(D0, D1);

    [HX_new, HZ_new] = coneAssembleOneZ(HX, HZ, zi, F1, F0, D1, D0);
    checkCSSCommutationZ(HX_new, HZ_new, 'coned HX_new,HZ_new');

    data.support = support;
    data.zSigns = zSigns;
    data.F1 = F1;
    data.F0 = F0;
    data.D1 = D1;
    data.D0 = D0;
    data.D1_pre = D1_pre;
    data.D0_pre = D0_pre;
    data.F0_pre = F0_pre;
    data.edges = edges;
    data.edges_pre = edges0;
    data.cellInfo = cellInfo;
    data.numVertices_Q = numel(support);
    data.numEdges_X_pre = size(D1_pre,1);
    data.numCycles_R_pre = size(D0_pre,1);
    data.numEdges_X = size(D1,1);
    data.numCycles_R = size(D0,1);
    data.numCellulationEdges = cellInfo.numAddedEdges;
    data.oldZWeight = nnz(HZ(zi,:));
    data.newZWeights = sum(HZ_new ~= 0, 2);
    data.newXWeights = sum(HX_new ~= 0, 2);
end
