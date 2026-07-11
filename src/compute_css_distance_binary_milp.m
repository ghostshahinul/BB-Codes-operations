function D = compute_css_distance_binary_milp(HX, HZ, opts)
%COMPUTE_CSS_DISTANCE_BINARY_MILP  Exact/time-limited CSS distance over F2.
%
% D = compute_css_distance_binary_milp(HX,HZ)
% D = compute_css_distance_binary_milp(HX,HZ,opts)
%
% INPUT
%   HX, HZ : CSS check matrices. Integer signed entries are allowed; the
%            function reduces them modulo 2 before computing distance.
%
% OPTIONS
%   opts.maxTimePerSolve   seconds for each MILP; default 300
%   opts.display           'off','iter','final'; default 'off'
%   opts.saveWitnesses     true/false; default true
%
% OUTPUT
%   D.dX, D.dZ             best values found
%   D.dFinal               min(D.dX,D.dZ)
%   D.dXExact, D.dZExact   true only if every required MILP was certified
%   D.exact                D.dXExact && D.dZExact
%   D.Xwitness, D.Zwitness minimum-weight candidates found
%   D.XlogicalBasis        basis of ker(HZ)/row(HX)
%   D.ZlogicalBasis        basis of ker(HX)/row(HZ)
%   D.Xsolves, D.Zsolves   per-MILP status
%
% MATHEMATICS
%   dZ = min wt(z), HX*z'=0, z not in row(HZ).
%   A nontrivial Z logical has odd overlap with at least one X-logical
%   basis vector. Therefore dZ is obtained from one MILP per X logical.
%
%   dX = min wt(x), HZ*x'=0, x not in row(HX).
%   A nontrivial X logical has odd overlap with at least one Z-logical
%   basis vector. Therefore dX is obtained from one MILP per Z logical.
%
% REQUIREMENT
%   MATLAB Optimization Toolbox (intlinprog).
%
% IMPORTANT
%   On a large qLDPC code, time-limited runs may produce upper bounds only.
%   D.exact is true only when all constituent optimization problems finish
%   with an optimality certificate.

    arguments
        HX
        HZ
        opts.maxTimePerSolve (1,1) double {mustBePositive} = 300
        opts.display (1,1) string = "off"
        opts.saveWitnesses (1,1) logical = true
    end

    if exist('intlinprog','file') ~= 2
        error(['intlinprog is unavailable. Install/enable MATLAB ' ...
               'Optimization Toolbox.']);
    end

    HX2 = sparse(mod(double(HX),2));
    HZ2 = sparse(mod(double(HZ),2));

    if size(HX2,2) ~= size(HZ2,2)
        error('HX and HZ must have the same number of columns.');
    end
    if nnz(mod(HX2*HZ2.',2)) ~= 0
        error('The matrices do not satisfy HX*HZ'' = 0 over F2.');
    end

    n = size(HX2,2);
    rankHX = gf2Rank(HX2);
    rankHZ = gf2Rank(HZ2);
    k = n-rankHX-rankHZ;

    fprintf('Binary CSS dimensions: n=%d, rank(HX)=%d, rank(HZ)=%d, k=%d\n', ...
            n,rankHX,rankHZ,k);

    fprintf('Constructing X-logical basis ker(HZ)/row(HX)...\n');
    LX = quotientLogicalBasis(HZ2,HX2);

    fprintf('Constructing Z-logical basis ker(HX)/row(HZ)...\n');
    LZ = quotientLogicalBasis(HX2,HZ2);

    if size(LX,1) ~= k || size(LZ,1) ~= k
        warning(['Logical basis size differs from k: size(LX,1)=%d, ' ...
                 'size(LZ,1)=%d, k=%d.'],size(LX,1),size(LZ,1),k);
    end

    fprintf('Solving %d MILPs for dZ using the X-logical basis...\n',size(LX,1));
    [dZ,Zwitness,Zinfo] = minLogicalMILPs(HX2,LX,opts);

    fprintf('Solving %d MILPs for dX using the Z-logical basis...\n',size(LZ,1));
    [dX,Xwitness,Xinfo] = minLogicalMILPs(HZ2,LZ,opts);

    D = struct();
    D.n = n;
    D.k = k;
    D.rankHX = rankHX;
    D.rankHZ = rankHZ;

    D.dX = dX;
    D.dZ = dZ;
    D.dFinal = min(dX,dZ);

    D.dXExact = Xinfo.allOptimal;
    D.dZExact = Zinfo.allOptimal;
    D.exact = D.dXExact && D.dZExact;

    D.XlogicalBasis = LX;
    D.ZlogicalBasis = LZ;
    D.Xsolves = Xinfo;
    D.Zsolves = Zinfo;

    if opts.saveWitnesses
        D.Xwitness = Xwitness;
        D.Zwitness = Zwitness;
    end

    fprintf('\nBest values found: dX=%g, dZ=%g, d_final=%g\n', ...
            D.dX,D.dZ,D.dFinal);

    if D.exact
        fprintf('All MILPs were certified optimal: d_final is exact.\n');
    else
        fprintf(['At least one MILP was not certified optimal. The finite ' ...
                 'values above are upper bounds, not a proven exact distance.\n']);
    end
end

function [best,bestWitness,info] = minLogicalMILPs(Hcomm,Lopposite,opts)
% Minimize a candidate logical that commutes with Hcomm and has odd overlap
% with one opposite-type logical basis vector.

    H = sparse(mod(double(Hcomm),2));
    L = sparse(mod(double(Lopposite),2));

    n = size(H,2);
    m = size(H,1);
    k = size(L,1);

    objectives = Inf(k,1);
    exitflags = NaN(k,1);
    relativeGaps = NaN(k,1);
    absoluteGaps = NaN(k,1);
    solveTimes = NaN(k,1);
    witnesses = cell(k,1);

    if k == 0
        best = Inf;
        bestWitness = sparse(1,n);
        info = struct('allOptimal',true,'objectives',objectives, ...
                      'exitflags',exitflags,'relativeGaps',relativeGaps, ...
                      'absoluteGaps',absoluteGaps,'solveTimes',solveTimes);
        return;
    end

    options = optimoptions('intlinprog', ...
        'Display',char(opts.display), ...
        'MaxTime',opts.maxTimePerSolve, ...
        'RelativeGapTolerance',0, ...
        'AbsoluteGapTolerance',0);

    rowWeights = full(sum(H~=0,2));
    tUpper = floor(rowWeights/2);

    for i = 1:k
        ell = sparse(double(L(i,:)));
        ellWeight = nnz(ell);
        if ellWeight == 0
            error('Encountered a zero logical-basis row.');
        end

        % Variables y=[v;t;s]:
        %   v in {0,1}^n is the candidate logical operator.
        %   H*v - 2*t = 0 enforces H*v=0 mod 2.
        %   ell*v - 2*s = 1 enforces odd logical overlap.
        Nvar = n+m+1;
        f = [ones(n,1); zeros(m+1,1)];
        intcon = 1:Nvar;

        Aeq = [H, -2*speye(m), sparse(m,1); ...
               ell, sparse(1,m), -2];
        beq = [zeros(m,1); 1];

        lb = zeros(Nvar,1);
        sUpper = floor((ellWeight-1)/2);
        ub = [ones(n,1); tUpper; sUpper];

        tic;
        [y,fval,exitflag,output] = intlinprog( ...
            f,intcon,[],[],Aeq,beq,lb,ub,options);
        solveTimes(i) = toc;

        exitflags(i) = exitflag;
        if isfield(output,'relativegap')
            relativeGaps(i) = output.relativegap;
        end
        if isfield(output,'absolutegap')
            absoluteGaps(i) = output.absolutegap;
        end

        if ~isempty(y)
            v = round(y(1:n));
            objectives(i) = round(fval);
            witnesses{i} = sparse(v(:).');
        end

        fprintf(['  logical class %2d/%2d: best=%g, exitflag=%g, ' ...
                 'time=%.1f s\n'], ...
                i,k,objectives(i),exitflags(i),solveTimes(i));
    end

    [best,idx] = min(objectives);
    if isfinite(best) && ~isempty(witnesses{idx})
        bestWitness = witnesses{idx};
    else
        bestWitness = sparse(1,n);
    end

    % With zero absolute and relative gap tolerances, exitflag==1 is the
    % conservative criterion used here for a certified optimum.
    allOptimal = all(exitflags == 1);

    info = struct();
    info.allOptimal = allOptimal;
    info.objectives = objectives;
    info.exitflags = exitflags;
    info.relativeGaps = relativeGaps;
    info.absoluteGaps = absoluteGaps;
    info.solveTimes = solveTimes;
end

function L = quotientLogicalBasis(Hcomm,Hstab)
% Basis representatives for ker(Hcomm)/row(Hstab) over F2.

    N = gf2Nullspace(Hcomm);
    [basisRows,pivots] = initializeGF2Basis(Hstab);
    logicals = false(0,size(Hcomm,2));

    for i = 1:size(N,1)
        original = full(logical(N(i,:)));
        [isNew,basisRows,pivots,reduced] = ...
            addToGF2Basis(original,basisRows,pivots);
        if isNew
            % reduced differs from the original by stabilizers and is
            % therefore a valid representative of the same logical class.
            logicals(end+1,:) = reduced; %#ok<AGROW>
        end
    end

    L = sparse(logicals);
end

function r = gf2Rank(A)
    [~,pivots] = initializeGF2Basis(A);
    r = numel(pivots);
end

function [basisRows,pivots] = initializeGF2Basis(A)
    basisRows = cell(0,1);
    pivots = zeros(0,1);

    for i = 1:size(A,1)
        [~,basisRows,pivots] = ...
            addToGF2Basis(full(logical(mod(A(i,:),2))),basisRows,pivots);
    end
end

function [isNew,basisRows,pivots,reduced] = ...
    addToGF2Basis(v,basisRows,pivots)

    v = logical(v);

    for j = 1:numel(pivots)
        if v(pivots(j))
            v = xor(v,basisRows{j});
        end
    end

    reduced = v;
    if ~any(v)
        isNew = false;
        return;
    end

    p = find(v,1,'first');
    insertAt = find(pivots>p,1,'first');

    if isempty(insertAt)
        pivots(end+1,1) = p;
        basisRows{end+1,1} = v;
    else
        pivots = [pivots(1:insertAt-1); p; pivots(insertAt:end)];
        basisRows = [basisRows(1:insertAt-1); {v}; ...
                     basisRows(insertAt:end)];
    end

    isNew = true;
end

function N = gf2Nullspace(A)
% Return a row basis of ker(A) over F2.

    A = logical(full(mod(A,2)));
    [m,n] = size(A);

    row = 1;
    pivots = zeros(0,1);

    for col = 1:n
        if row > m
            break;
        end

        idx = find(A(row:end,col),1,'first');
        if isempty(idx)
            continue;
        end
        idx = idx+row-1;

        if idx ~= row
            tmp = A(row,:);
            A(row,:) = A(idx,:);
            A(idx,:) = tmp;
        end

        other = find(A(:,col));
        other(other==row) = [];
        for q = other.'
            A(q,:) = xor(A(q,:),A(row,:));
        end

        pivots(end+1,1) = col; %#ok<AGROW>
        row = row+1;
    end

    freeColumns = setdiff(1:n,pivots);
    N = false(numel(freeColumns),n);

    for j = 1:numel(freeColumns)
        f = freeColumns(j);
        N(j,f) = true;

        for r = 1:numel(pivots)
            if A(r,f)
                N(j,pivots(r)) = true;
            end
        end
    end

    N = sparse(N);
end
