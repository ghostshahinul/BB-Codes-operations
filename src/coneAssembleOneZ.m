function [HX_new, HZ_new] = coneAssembleOneZ(HX, HZ, zi, F1, F0, D1, D0)
% coneAssembleOneZ  Assemble one-row coned matrices over Z.
% Column order: [new edge qubits | old qubits]
    HX = double(HX);
    HZ = double(HZ);
    n = size(HX,2);
    E = size(D1,1);

    HZ_r = HZ;
    HZ_r(zi,:) = [];

    HX_new = [sparse(D0), sparse(size(D0,1), n); sparse(F0), sparse(HX)];
    HZ_new = [sparse(D1.'), sparse(F1.'); sparse(size(HZ_r,1), E), sparse(HZ_r)];
end
