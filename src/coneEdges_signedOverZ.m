function edges = coneEdges_signedOverZ(HX, HZ, zi)
% coneEdges_signedOverZ  Build Xi edge list from signed integer overlaps.
% For each X row r, b = HX(r,support).*HZ(zi,support), sum(b)=0.
% Pair positive and negative entries to create signed edges.

    HX = double(HX);
    HZ = double(HZ);
    [nX, n] = size(HX); %#ok<ASGLU>
    if size(HZ,2) ~= n
        error('HX and HZ must have the same number of columns.');
    end

    [support, ~] = coneSupportZ(HZ, zi);
    z = HZ(zi, support);

    owner = [];
    uLocal = [];
    vLocal = [];
    uGlobal = [];
    vGlobal = [];
    amount = [];

    for r = 1:size(HX,1)
        b = HX(r, support) .* z;
        if all(b == 0)
            continue;
        end
        if sum(b) ~= 0
            error('Row %d has nonzero integer overlap with Z row %d: sum = %g.', r, zi, sum(b));
        end

        pos = find(b > 0);
        neg = find(b < 0);
        posAmt = b(pos);
        negAmt = -b(neg);

        ip = 1;
        ineg = 1;
        while ip <= numel(pos) && ineg <= numel(neg)
            a = min(posAmt(ip), negAmt(ineg));

            owner(end+1,1) = r;
            uLocal(end+1,1) = pos(ip);      % D1 gets -a here
            vLocal(end+1,1) = neg(ineg);    % D1 gets +a here
            uGlobal(end+1,1) = support(pos(ip));
            vGlobal(end+1,1) = support(neg(ineg));
            amount(end+1,1) = a;

            posAmt(ip) = posAmt(ip) - a;
            negAmt(ineg) = negAmt(ineg) - a;

            if posAmt(ip) == 0
                ip = ip + 1;
            end
            if negAmt(ineg) == 0
                ineg = ineg + 1;
            end
        end

        if ip <= numel(pos) || ineg <= numel(neg)
            error('Internal flow decomposition failed on X row %d.', r);
        end
    end

    edges.ownerXRow = owner;
    edges.uLocal = uLocal;
    edges.vLocal = vLocal;
    edges.uGlobal = uGlobal;
    edges.vGlobal = vGlobal;
    edges.amount = amount;
    edges.isCellulation = false(numel(owner),1);
    edges.numEdges = numel(owner);
    edges.support = support;
end
