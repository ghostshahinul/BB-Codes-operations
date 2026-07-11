function stats = coneStatsZ(HX, HZ)
% coneStatsZ  Basic row/column weights and entry range.
    HX = double(HX);
    HZ = double(HZ);
    stats.sizeHX = size(HX);
    stats.sizeHZ = size(HZ);
    stats.wX = max(sum(HX ~= 0, 2));
    stats.wZ = max(sum(HZ ~= 0, 2));
    stats.qX = max(sum(HX ~= 0, 1));
    stats.qZ = max(sum(HZ ~= 0, 1));
    stats.uniqueHX = unique(HX(:));
    stats.uniqueHZ = unique(HZ(:));
    stats.nnzCommutation = nnz(HX*HZ.');
    stats.maxAbsCommutation = max(abs(HX*HZ.'), [], 'all');
end
