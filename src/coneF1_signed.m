function F1 = coneF1_signed(HZ, zi)
% coneF1_signed  Signed embedding f1: Q_i -> old qubits C1.
    HZ = double(HZ);
    [support, zSigns] = coneSupportZ(HZ, zi);
    n = size(HZ,2);
    w = numel(support);
    F1 = sparse(support, 1:w, zSigns, n, w);
end
