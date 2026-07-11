function r = ufFind(parent, x)
% ufFind  Find root in union-find parent array.
    while parent(x) ~= x
        x = parent(x);
    end
    r = x;
end
