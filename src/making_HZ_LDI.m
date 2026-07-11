function [HZ_ldi, info] = make_HZ_LDI(HX, HZ, varargin)
%MAKE_HZ_LDI  Create a signed Z-check matrix HZ_ldi with integer commutation.
%
%   [HZ_ldi, info] = make_HZ_LDI(HX, HZ)
%
%   Given binary CSS check matrices HX and HZ with the same number of
%   columns, this function tries to replace every nonzero entry of HZ by
%   +1 or -1 so that
%
%       HX * HZ_ldi.' == 0
%
%   over the integers, while preserving the binary support:
%
%       mod(HZ_ldi,2) == HZ.
%
%   This is the signed-CSS LDI seed construction used for BB-code coning.
%
%   IMPORTANT:
%   - HZ alone is not enough. The signing depends on overlaps with HX.
%   - The fast graph method works when every nonzero X-Z overlap has size 2.
%   - For larger even overlaps, the function uses backtracking per Z row.
%   - If no +/- signing exists for some row, the function throws an error.
%
%   Optional name-value arguments:
%       'MaxBacktrackWeight' : maximum Z-row weight for backtracking.
%                              Default = 28.
%       'Verbose'            : print row progress. Default = false.
%
%   Example:
%       HX = readmatrix("BB72_HX.csv");
%       HZ = readmatrix("BB72_HZ.csv");
%       [HZ_ldi, info] = make_HZ_LDI(HX, HZ);
%       nnz(double(HX)*double(HZ_ldi).')
%       unique(HZ_ldi(:))
%       isequal(mod(HZ_ldi,2), double(HZ~=0))
%
%   Author: generated for BB-code / LDI calculations.

    p = inputParser;
    addParameter(p, 'MaxBacktrackWeight', 28, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Verbose', false, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    maxBacktrackWeight = p.Results.MaxBacktrackWeight;
    verbose = logical(p.Results.Verbose);

    HX = double(HX ~= 0);
    HZ_bin = double(HZ ~= 0);

    [mX, nXcols] = size(HX);
    [mZ, nZcols] = size(HZ_bin);

    if nXcols ~= nZcols
        error('HX and HZ must have the same number of columns. Got %d and %d.', nXcols, nZcols);
    end

    HZ_ldi = zeros(mZ, nZcols);

    info = struct();
    info.success = false;
    info.methods = strings(mZ,1);
    info.supports = cell(mZ,1);
    info.signs = cell(mZ,1);
    info.constraintWeights = cell(mZ,1);

    for j = 1:mZ
        supp = find(HZ_bin(j,:) ~= 0);
        w = numel(supp);
        info.supports{j} = supp;

        if verbose
            fprintf('Signing Z row %d/%d, weight %d...\n', j, mZ, w);
        end

        if w == 0
            info.methods(j) = "empty";
            info.signs{j} = [];
            continue;
        end

        % Build constraint matrix A. Each row says:
        %       sum_{c in overlap} sign(c) = 0.
        % Here columns of A are local support coordinates inside supp.
        A = [];
        for i = 1:mX
            loc = find(HX(i, supp) ~= 0);
            r = numel(loc);

            if r == 0
                continue;
            end

            if mod(r,2) == 1
                error(['No +/- LDI signing possible for Z row %d: ', ...
                       'X row %d overlaps it in odd size %d.'], j, i, r);
            end

            row = zeros(1,w);
            row(loc) = 1;
            A = [A; row]; %#ok<AGROW>
        end

        info.constraintWeights{j} = sum(A ~= 0, 2);

        if isempty(A)
            s = ones(1,w);
            info.methods(j) = "unconstrained";
        elseif all(sum(A ~= 0, 2) == 2)
            [s, ok] = solveByBipartiteGraph(A);
            if ~ok
                error('No +/- LDI signing possible for Z row %d: conflict graph is not bipartite.', j);
            end
            info.methods(j) = "bipartite_graph";
        else
            if w > maxBacktrackWeight
                error(['Z row %d has weight %d, which exceeds MaxBacktrackWeight=%d. ', ...
                       'Increase MaxBacktrackWeight or use an ILP/SAT solver.'], j, w, maxBacktrackWeight);
            end

            [s, ok] = solveByBacktracking(A);
            if ~ok
                error('No +/- LDI signing found for Z row %d by backtracking.', j);
            end
            info.methods(j) = "backtracking";
        end

        HZ_ldi(j, supp) = s;
        info.signs{j} = s;
    end

    % Final verification.
    C = double(HX) * double(HZ_ldi).';
    info.checkMatrix = C;
    info.nnzCheck = nnz(C);
    info.maxAbsCheck = max(abs(C(:)));
    info.sameMod2Support = isequal(mod(HZ_ldi,2), HZ_bin);
    info.uniqueEntries = unique(HZ_ldi(:));
    info.success = (info.nnzCheck == 0) && info.sameMod2Support;

    if ~info.success
        error('Internal verification failed: HX*HZ_ldi.'' is not zero or mod-2 support changed.');
    end
end

function [s, ok] = solveByBipartiteGraph(A)
% A has one row per edge, two ones per row.
    w = size(A,2);
    Adj = false(w,w);

    for e = 1:size(A,1)
        loc = find(A(e,:) ~= 0);
        a = loc(1);
        b = loc(2);
        Adj(a,b) = true;
        Adj(b,a) = true;
    end

    color = zeros(1,w);

    for start = 1:w
        if color(start) ~= 0
            continue;
        end

        color(start) = 1;
        queue = start;

        while ~isempty(queue)
            u = queue(1);
            queue(1) = [];
            nbrs = find(Adj(u,:));

            for v = nbrs
                if color(v) == 0
                    color(v) = -color(u);
                    queue(end+1) = v; %#ok<AGROW>
                elseif color(v) == color(u)
                    s = [];
                    ok = false;
                    return;
                end
            end
        end
    end

    s = color;
    ok = all(A*s.' == 0);
end

function [s, ok] = solveByBacktracking(A)
% Search signs s in {+1,-1}^w satisfying A*s.' = 0.
% Uses basic feasibility pruning.
    w = size(A,2);
    deg = sum(A ~= 0, 1);
    [~, order] = sort(deg, 'descend');

    s = zeros(1,w);
    ok = recurse(1);

    function tf = recurse(t)
        if t > w
            tf = all(A*s.' == 0);
            return;
        end

        v = order(t);

        for val = [1, -1]
            s(v) = val;
            if feasible()
                if recurse(t+1)
                    tf = true;
                    return;
                end
            end
        end

        s(v) = 0;
        tf = false;
    end

    function tf = feasible()
        partial = A*s.';
        remaining = A*(s.' == 0);

        % Can the remaining +/- variables cancel the current partial sum?
        tf = all(abs(partial) <= remaining) && ...
             all(mod(remaining - abs(partial), 2) == 0);
    end
end

HX = [
1 0 0 1 1 0 0 1 0 0 0 0;
0 0 0 0 0 0 0 1 0 1 0 0;
0 0 0 1 0 0 1 0 0 0 0 0;
0 1 0 0 0 0 1 0 1 1 0 0;
1 1 0 0 0 0 0 0 0 0 0 1;
0 0 0 0 0 1 0 0 0 0 0 1;
0 0 0 0 0 0 0 0 1 0 1 0;
0 0 1 0 0 0 0 0 0 0 1 0;
1 0 1 0 0 0 0 0 0 0 0 0
];

HZ = [
0 1 0 1 1 1 1 0 0 0 0 1;
1 0 1 0 0 1 0 1 1 1 1 1;
0 0 0 1 0 0 1 1 0 1 0 0
];


[HZ_ldi, info] = make_HZ_LDI(HX, HZ);
C = double(HX) * double(HZ_ldi).';

nnz(C)
max(abs(C(:)))
unique(HZ_ldi(:))
isequal(mod(HZ_ldi,2), double(HZ ~= 0))