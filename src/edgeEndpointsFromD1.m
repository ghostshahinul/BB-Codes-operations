function [u, v] = edgeEndpointsFromD1(D1, e)
% edgeEndpointsFromD1  Return the two endpoint vertices of edge row e.
    idx = find(D1(e,:) ~= 0);
    if numel(idx) ~= 2
        error('D1 row %d must have exactly two nonzero entries.', e);
    end
    u = idx(1);
    v = idx(2);
end
