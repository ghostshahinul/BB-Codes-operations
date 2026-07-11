function [tail, head, amount] = edgeOrientationFromD1(D1, e)
% edgeOrientationFromD1
% D1 row should have one negative and one positive entry.
% tail is the negative endpoint, head is the positive endpoint.

    idx = find(D1(e,:) ~= 0);
    if numel(idx) ~= 2
        error('D1 row %d must have exactly two nonzero entries.', e);
    end

    vals = D1(e,idx);
    neg = find(vals < 0);
    pos = find(vals > 0);

    if numel(neg) ~= 1 || numel(pos) ~= 1
        error('D1 row %d must have one negative and one positive endpoint.', e);
    end

    amountNeg = abs(vals(neg));
    amountPos = abs(vals(pos));
    if amountNeg ~= amountPos
        error('D1 row %d has unmatched endpoint magnitudes.', e);
    end

    tail = idx(neg);
    head = idx(pos);
    amount = amountNeg;
end
