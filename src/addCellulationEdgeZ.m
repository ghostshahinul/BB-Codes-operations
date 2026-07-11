function [D1_new, F0_new, edges_new, newEdgeIdx] = addCellulationEdgeZ(D1_old, F0_old, edges_old, tailLocal, headLocal)
% addCellulationEdgeZ
% Append one auxiliary edge tailLocal -> headLocal to D1.
% It has D1(new,tailLocal)=-1 and D1(new,headLocal)=+1.
% Its F0 column is zero because it does not come from an old X-check.

    D1_new = D1_old;
    Eold = size(D1_old,1);
    V = size(D1_old,2);
    newEdgeIdx = Eold + 1;

    newRow = zeros(1,V);
    newRow(tailLocal) = -1;
    newRow(headLocal) = +1;
    D1_new(newEdgeIdx,:) = newRow;

    F0_new = F0_old;
    F0_new(:,newEdgeIdx) = 0;

    edges_new = edges_old;
    edges_new.ownerXRow(newEdgeIdx,1) = 0;
    edges_new.uLocal(newEdgeIdx,1) = tailLocal;
    edges_new.vLocal(newEdgeIdx,1) = headLocal;

    if isfield(edges_new, 'support') && ~isempty(edges_new.support)
        edges_new.uGlobal(newEdgeIdx,1) = edges_new.support(tailLocal);
        edges_new.vGlobal(newEdgeIdx,1) = edges_new.support(headLocal);
    else
        edges_new.uGlobal(newEdgeIdx,1) = NaN;
        edges_new.vGlobal(newEdgeIdx,1) = NaN;
    end

    edges_new.amount(newEdgeIdx,1) = 1;
    edges_new.isCellulation(newEdgeIdx,1) = true;
    edges_new.numEdges = newEdgeIdx;
end
