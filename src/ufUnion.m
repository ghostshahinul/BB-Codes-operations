function parent = ufUnion(parent, a, b)
% ufUnion  Union roots a and b.
    ra = ufFind(parent, a);
    rb = ufFind(parent, b);
    if ra ~= rb
        parent(rb) = ra;
    end
end
