function [support, zSigns] = coneSupportZ(HZ, zi)
% coneSupportZ  Support and signed entries of selected Z row.
    HZ = double(HZ);
    if zi < 1 || zi > size(HZ,1)
        error('zi is outside the row range of HZ.');
    end
    support = find(HZ(zi,:) ~= 0);
    zSigns = HZ(zi, support);
    if isempty(support)
        error('Selected Z row has empty support.');
    end
end
