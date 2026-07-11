function R = run_BB72_conditional_weight_reduction_Z(opts)
%RUN_BB72_CONDITIONAL_WEIGHT_REDUCTION_Z
% Conditional Hastings-style pipeline for the BB [[72,12,6]] code:
%   1) build BB matrices
%   2) construct signed HZ_LDI over Z
%   3) copying only if qX>3
%   4) gauging only if wX>3
%   5) thickening/heights only if qZ>3 or explicitly requested
%   6) simultaneous coning + cellulation over Z
%   7) report/save HX, HZ, n, nX, nZ, wX, wZ, qX, qZ
%   8) optionally attempt exact binary distance using MILP
%
% IMPORTANT
% ---------
% This code performs exact INTEGER arithmetic.  It never uses mod 2 while
% constructing the matrices.  At every stage it verifies
%
%                  HX * HZ.' == 0 over Z.
%
% The paper is written over F2.  The signs below implement the required
% nonbinary/integer sign conventions:
%   - copying repetition checks use [ +1  -1 ]
%   - gauging auxiliary columns use oriented +1/-1 incidence
%   - thickening uses the signed tensor-product differential
%   - the mapping-cone minus is absorbed into D1 exactly once
%
% INPUT
% -----
% opts is optional.  Useful fields:
%
%   opts.targetQX/targetWX/targetQZ default 3
%   opts.maxCycleLen          maximum D0 row weight after cellulation,
%                             default 4
%   opts.coneOnlyRowsAbove    cone selected-height Z rows whose weight is
%                             above this value, default 5
%   opts.writeCSV             write dense CSV files, default false
%   opts.saveMAT              save MAT result, default true
%   opts.outputPrefix         default 'BB72_full_WR_Z'
%   opts.distanceMode         'report_only' (default) or 'milp_binary'
%   opts.distanceTimePerSolve seconds per MILP, default 60
%   opts.verbose              default true
%
% OUTPUT
% ------
% R.HX_final, R.HZ_final      final sparse integer matrices
% R.final                     final statistics
% R.stages                    statistics at every stage
% R.distance                  distance report / optional MILP result
% R.heights                   selected heights
% R.copy, R.gauge, R.thicken, R.cone  metadata
%
% EXAMPLE
% -------
%   opts = struct;
%   opts.targetQX = 3;
%   opts.targetWX = 3;
%   opts.targetQZ = 3;
%   opts.maxCycleLen = 4;
%   opts.distanceMode = 'report_only';
%   R = run_BB72_conditional_weight_reduction_Z(opts);
%
% To attempt an exact qubit distance calculation (potentially very slow):
%   opts.distanceMode = 'milp_binary';
%   opts.distanceTimePerSolve = 120;
%   R = run_BB72_conditional_weight_reduction_Z(opts);
%
% The exact distance of the large final code is a separate optimization
% problem.  The matrix construction itself cannot infer it from weights.

    if nargin < 1 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    if opts.verbose
        fprintf('\n============================================================\n');
        fprintf('BB [[72,12,6]] CONDITIONAL weight reduction over Z\n');
        fprintf('targets: qX<=%d, wX<=%d, qZ<=%d; maxCycleLen=%d\n', ...
            opts.targetQX,opts.targetWX,opts.targetQZ,opts.maxCycleLen);
        fprintf('============================================================\n');
    end

    %% 0. Build BB [[72,12,6]] matrices and signed HZ_LDI
    [HX0, HZ_binary, bbMeta] = buildBB72();
    HZ0 = signHZ_LDI_bipartite(HX0, HZ_binary);
    assertCommutesZ(HX0, HZ0, 'signed BB seed');

    R = struct();
    R.options = opts;
    R.bb = bbMeta;
    R.HX_seed = HX0;
    R.HZ_binary = HZ_binary;
    R.HZ_LDI_seed = HZ0;
    R.stages.seed = codeStatsZ(HX0, HZ0, 'signed seed');
    printStats(R.stages.seed, opts.verbose);

    %% 1. Copying -- ONLY when qX exceeds the target
    s0 = R.stages.seed;
    needCopying = opts.forceCopying || (s0.qX > opts.targetQX);

    if needCopying
        [HX1, HZ1, copyMeta] = copyingOverZ(HX0, HZ0);
        copyMeta.applied = true;
        copyMeta.reason = sprintf('qX=%d exceeds target qX=%d.',s0.qX,opts.targetQX);
        stageLabel = 'after copying';
    else
        HX1 = HX0;
        HZ1 = HZ0;
        copyMeta = struct();
        copyMeta.applied = false;
        copyMeta.qXcopies = 1;
        copyMeta.reason = sprintf('Skipped: qX=%d is already <= %d.',s0.qX,opts.targetQX);
        stageLabel = 'copying skipped';
    end
    assertCommutesZ(HX1, HZ1, stageLabel);
    R.copy = copyMeta;
    R.stages.copying = codeStatsZ(HX1, HZ1, stageLabel);
    printDecision('Copying',copyMeta.applied,copyMeta.reason,opts.verbose);
    printStats(R.stages.copying, opts.verbose);

    %% 2. Gauging -- ONLY rows whose weight exceeds targetWX
    s1 = R.stages.copying;
    needGauging = opts.forceGauging || (s1.wX > opts.targetWX);

    if needGauging
        [HX2, HZ2, gaugeMeta] = gaugeAllXOverZ(HX1, HZ1, opts.targetWX);
        gaugeMeta.applied = true;
        gaugeMeta.reason = sprintf('wX=%d exceeds target wX=%d.',s1.wX,opts.targetWX);
        stageLabel = 'after gauging';
    else
        HX2 = HX1;
        HZ2 = HZ1;
        gaugeMeta = struct();
        gaugeMeta.applied = false;
        gaugeMeta.numGaugeQubits = 0;
        gaugeMeta.reason = sprintf('Skipped: wX=%d is already <= %d.',s1.wX,opts.targetWX);
        stageLabel = 'gauging skipped';
    end
    assertCommutesZ(HX2, HZ2, stageLabel);
    R.gauge = gaugeMeta;
    R.stages.gauging = codeStatsZ(HX2, HZ2, stageLabel);
    printDecision('Gauging',gaugeMeta.applied,gaugeMeta.reason,opts.verbose);
    printStats(R.stages.gauging, opts.verbose);

    useCertifiedPlan = false;

    %% 3. Thickening/choosing heights -- only when qZ needs reduction
    % Thickening may also be requested explicitly for distance balancing.
    % Exact dX is not available from row/column weights alone, so the code
    % cannot honestly trigger on "dX decreased" without a distance solver.
    s2 = R.stages.gauging;
    needThickening = opts.forceThickening || opts.balanceAfterGauging || ...
        (s2.qZ > opts.targetQZ);

    if needThickening
        % For this exact conditional BB route, ell=7 and the following
        % height assignment are certified to achieve qZ=3.  If the input
        % changes, the generic MILP search is used instead.
        [useCertifiedPlan, certifiedEll, certifiedHeights] = ...
            certifiedConditionalBB72Heights(HX2,HZ2,opts.targetQZ);

        if opts.useCertifiedBB72Plan && useCertifiedPlan
            ellUsed = certifiedEll;
            fixedHeights = certifiedHeights;
            heightSearch = struct('method','certified BB72 conditional plan', ...
                'testedEll',2:ellUsed,'feasible',true);
        else
            [ellUsed,fixedHeights,heightSearch] = findMinimalHeightPlanOverZ( ...
                HX2,HZ2,opts.targetQZ,opts.maxEll,opts.heightMILPTime);
        end

        [HX3, HZ3, thickMeta] = thickenChooseHeightsOverZ( ...
            HX2,HZ2,ellUsed,fixedHeights);
        thickMeta.applied = true;
        thickMeta.reason = sprintf(['Applied because qZ=%d exceeds target %d ' ...
            'or distance balancing was requested.'],s2.qZ,opts.targetQZ);
        thickMeta.heightSearch = heightSearch;
        stageLabel = 'after thickening/heights';
    else
        HX3 = HX2;
        HZ3 = HZ2;
        ellUsed = 1;
        thickMeta = struct();
        thickMeta.applied = false;
        thickMeta.ell = 1;
        thickMeta.heights = [];
        thickMeta.selectedTopRowIndices = (1:size(HZ3,1)).';
        thickMeta.bottomRowIndices = zeros(0,1);
        thickMeta.reason = sprintf('Skipped: qZ=%d is already <= %d.',s2.qZ,opts.targetQZ);
        stageLabel = 'thickening skipped';
    end
    assertCommutesZ(HX3, HZ3, stageLabel);
    R.thicken = thickMeta;
    R.heights = thickMeta.heights;
    R.ellUsed = ellUsed;
    R.stages.thickening = codeStatsZ(HX3, HZ3, stageLabel);
    printDecision('Thickening/heights',thickMeta.applied,thickMeta.reason,opts.verbose);
    printStats(R.stages.thickening, opts.verbose);

    if thickMeta.applied && R.stages.thickening.qZ > opts.targetQZ
        error('Height selection failed: obtained qZ=%d, target was %d.', ...
            R.stages.thickening.qZ,opts.targetQZ);
    end

    %% 4. Coning -- ONLY high-weight Z rows
    rowWeightsZ3 = full(sum(HZ3 ~= 0,2));
    if thickMeta.applied
        candidateRows = thickMeta.selectedTopRowIndices;
    else
        candidateRows = (1:size(HZ3,1)).';
    end
    coneRows = candidateRows(rowWeightsZ3(candidateRows) > opts.coneOnlyRowsAbove);

    if isempty(coneRows)
        HX4 = HX3;
        HZ4 = HZ3;
        coneMeta = struct([]);
        coneApplied = false;
        coneReason = sprintf('Skipped: no candidate Z row has weight > %d.', ...
            opts.coneOnlyRowsAbove);
    else
        [HX4, HZ4, coneMeta] = coneManyZOverZCellulated( ...
            HX3,HZ3,coneRows,opts.maxCycleLen, ...
            opts.targetCellulatedVertexDegree,opts.maxCellulationTries,opts.verbose);
        coneApplied = true;
        coneReason = sprintf('Coned %d Z rows with weight > %d.', ...
            numel(coneRows),opts.coneOnlyRowsAbove);
    end
    assertCommutesZ(HX4, HZ4, 'final coned code');
    printDecision('Coning',coneApplied,coneReason,opts.verbose);

    R.coneRows = coneRows;
    R.cone = coneMeta;
    R.HX_final = HX4;
    R.HZ_final = HZ4;
    R.final = codeStatsZ(HX4, HZ4, 'final');
    R.stages.final = R.final;
    printStats(R.final, opts.verbose);

    % Certified result for the current BB [[72,12,6]] conditional route.
    if opts.verifyConditionalBB72Target && useCertifiedPlan
        target = [3938,1886,2052,6,6,5,3];
        actual = [R.final.n,R.final.nX,R.final.nZ, ...
            R.final.wX,R.final.wZ,R.final.qX,R.final.qZ];
        if ~isequal(actual,target)
            error(['Conditional BB72 verification failed. Expected ' ...
                '[n nX nZ wX wZ qX qZ]=[%s], obtained [%s].'], ...
                num2str(target),num2str(actual));
        end
    end

    %% Dimension checks over selected primes (optional but useful for LDI)
    if opts.computePrimeRanks
        primes = opts.primes;
        dimRows = zeros(numel(primes), 4);
        for a = 1:numel(primes)
            p = primes(a);
            rX = rankModPrime(HX4, p);
            rZ = rankModPrime(HZ4, p);
            k = size(HX4,2) - rX - rZ;
            dimRows(a,:) = [p, rX, rZ, k];
        end
        R.dimensionByPrime = array2table(dimRows, ...
            'VariableNames', {'p','rankHX','rankHZ','k'});
        if opts.verbose
            disp(R.dimensionByPrime);
        end
    else
        R.dimensionByPrime = [];
    end

    %% Distance report / optional exact binary MILP
    opts.ell = ellUsed;
    R.distance = distanceReport(HX4, HZ4, opts, copyMeta);

    %% Save outputs
    if opts.saveMAT
        matName = [opts.outputPrefix '_results.mat'];
        save(matName, 'R', '-v7.3');
        if opts.verbose
            fprintf('Saved MAT file: %s\n', matName);
        end
    end

    if opts.writeCSV
        hxName = [opts.outputPrefix '_HX_final.csv'];
        hzName = [opts.outputPrefix '_HZ_final.csv'];
        warning(['Writing dense CSV files for matrices of this size can be ' ...
                 'slow and consume substantial disk space.']);
        writematrix(full(HX4), hxName);
        writematrix(full(HZ4), hzName);
        if opts.verbose
            fprintf('Saved CSV files: %s, %s\n', hxName, hzName);
        end
    end

    if opts.verbose
        fprintf('\nFinal integer commutator: nnz(HX*HZ'') = %d\n', nnz(HX4*HZ4.'));
        fprintf('Finished.\n\n');
    end
end

%% ========================================================================
%  OPTIONS
%  ========================================================================
function opts = applyDefaults(opts)
    defaults.targetQX = 3;
    defaults.targetWX = 3;
    defaults.targetQZ = 3;
    defaults.coneOnlyRowsAbove = 5;
    defaults.maxCycleLen = 4;
    defaults.targetCellulatedVertexDegree = 5; % gives replacement wZ <= 6
    defaults.maxCellulationTries = 30000;
    defaults.maxEll = 14;
    defaults.heightMILPTime = 60;
    defaults.forceCopying = false;
    defaults.forceGauging = false;
    defaults.forceThickening = false;
    defaults.balanceAfterGauging = false;
    defaults.useCertifiedBB72Plan = true;
    defaults.verifyConditionalBB72Target = true;
    defaults.writeCSV = false;
    defaults.saveMAT = true;
    defaults.outputPrefix = 'BB72_conditional_WR_Z';
    defaults.distanceMode = 'report_only';
    defaults.distanceTimePerSolve = 60;
    defaults.verbose = true;
    defaults.computePrimeRanks = true;
    defaults.primes = [2 3 5 7];

    names = fieldnames(defaults);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(opts,name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end

    if opts.targetQX < 1 || opts.targetWX ~= 3 || opts.targetQZ < 1
        error('This implementation requires targetWX=3 and positive targets.');
    end
    if opts.maxCycleLen ~= 4
        error('Balanced cellulation currently supports maxCycleLen=4.');
    end
    if opts.targetCellulatedVertexDegree < 3
        error('targetCellulatedVertexDegree must be at least 3.');
    end
end

function printDecision(name,applied,reason,verbose)
    if ~verbose, return; end
    if applied, action='APPLIED'; else, action='SKIPPED'; end
    fprintf('%-27s %s -- %s\n',name,action,reason);
end

%% ========================================================================
%  BB [[72,12,6]] CONSTRUCTION
%  ========================================================================
function [HX, HZ, meta] = buildBB72()
% Uses l=m=6, A=x^3+y+y^2, B=y^3+x+x^2.
% The shift orientation is chosen to match the prior checkpoint:
% row 1 of HZ has zero-based support {3,24,30,40,41,54}.

    ell = 6;
    m = 6;
    N = ell*m;

    X = bbShiftMatrix(ell, m, 'x', 1);
    Y = bbShiftMatrix(ell, m, 'y', 1);

    A = X^3 + Y + Y^2;
    B = Y^3 + X + X^2;
    A = sparse(mod(A,2));
    B = sparse(mod(B,2));

    HX = sparse([A, B]);
    HZ = sparse([B.', A.']);

    meta.ell = ell;
    meta.m = m;
    meta.n = 2*N;
    meta.A = A;
    meta.B = B;
    meta.polynomialA = 'x^3 + y + y^2';
    meta.polynomialB = 'y^3 + x + x^2';
end

function P = bbShiftMatrix(ell, m, axisName, power)
% Column basis vector (a,b) is shifted by -power.  This convention matches
% the previously generated BB72 matrices and support indexing.

    N = ell*m;
    rows = zeros(N,1);
    cols = (1:N).';

    t = 1;
    for a = 0:ell-1
        for b = 0:m-1
            if strcmpi(axisName,'x')
                aa = mod(a-power, ell);
                bb = b;
            elseif strcmpi(axisName,'y')
                aa = a;
                bb = mod(b-power, m);
            else
                error('axisName must be x or y.');
            end
            rows(t) = aa*m + bb + 1;
            t = t + 1;
        end
    end

    P = sparse(rows, cols, 1, N, N);
end

function HZsigned = signHZ_LDI_bipartite(HX, HZbinary)
% Sign every binary Z row by bipartite coloring of its conflict graph.
% The output has entries in {-1,0,+1} and satisfies HX*HZsigned.'=0 over Z.

    HX = sparse(double(HX));
    HZbinary = sparse(double(HZbinary ~= 0));
    [nZ,n] = size(HZbinary);
    nX = size(HX,1);
    HZsigned = sparse(nZ,n);

    for j = 1:nZ
        supp = find(HZbinary(j,:) ~= 0);
        w = numel(supp);
        Adj = sparse(w,w);

        for r = 1:nX
            local = find(HX(r,supp) ~= 0);
            if numel(local) == 2
                Adj(local(1),local(2)) = 1;
                Adj(local(2),local(1)) = 1;
            elseif ~isempty(local)
                error(['Signing graph failed for Z row %d: X row %d has ' ...
                       'overlap size %d instead of 0 or 2.'], j, r, numel(local));
            end
        end

        color = zeros(1,w);
        for v = 1:w
            if color(v) ~= 0
                continue;
            end
            color(v) = 1;
            queue = v;
            while ~isempty(queue)
                u = queue(1);
                queue(1) = [];
                nbrs = find(Adj(u,:) ~= 0);
                for q = nbrs
                    if color(q) == 0
                        color(q) = -color(u);
                        queue(end+1) = q; %#ok<AGROW>
                    elseif color(q) == color(u)
                        error('Conflict graph for Z row %d is not bipartite.', j);
                    end
                end
            end
        end
        HZsigned(j,supp) = color;
    end

    assertCommutesZ(HX, HZsigned, 'signed HZ construction');
end

%% ========================================================================
%  COPYING OVER Z
%  ========================================================================
function [HXnew, HZnew, meta] = copyingOverZ(HX, HZ)
% Original Hastings copying with global qX copies per qubit.
% Over Z, the repetition-link rows are oriented [1,-1], not [1,1].

    HX = sparse(double(HX));
    HZ = sparse(double(HZ));
    [nX,n] = size(HX);
    qX = max(full(sum(HX ~= 0,1)));

    if qX < 1
        error('HX has no nonzero columns.');
    end

    HXnew = sparse(nX + n*(qX-1), n*qX);

    % Copy each old X occurrence into a different copy of that qubit.
    for col = 1:n
        rows = find(HX(:,col) ~= 0);
        vals = full(HX(rows,col));
        if numel(rows) > qX
            error('Internal copying error: column degree exceeds qX.');
        end
        for t = 1:numel(rows)
            newCol = (col-1)*qX + t;
            HXnew(rows(t),newCol) = vals(t);
        end
    end

    % Oriented repetition checks between adjacent copies.
    rowPtr = nX;
    for col = 1:n
        for c = 1:qX-1
            rowPtr = rowPtr + 1;
            c1 = (col-1)*qX + c;
            c2 = c1 + 1;
            HXnew(rowPtr,c1) = +1;
            HXnew(rowPtr,c2) = -1;
        end
    end

    % Duplicate each signed Z entry on every copy.
    HZnew = kron(HZ, sparse(ones(1,qX)));

    meta.qXcopies = qX;
    meta.oldXRows = 1:nX;
    meta.linkXRows = nX+1:size(HXnew,1);
    meta.copyColumnMap = reshape(1:n*qX, qX, n).';
end

%% ========================================================================
%  GAUGING OVER Z
%  ========================================================================
function [HXnew, HZnew, meta] = gaugeAllXOverZ(HX, HZ, targetWeight)
% Replace every X row of weight >3 by an oriented path gadget.
% The new Z entries are solved recursively over Z so that every new X row
% commutes exactly with every Z row.

    if targetWeight ~= 3
        error('This implementation currently supports targetWeight=3 only.');
    end

    HX = sparse(double(HX));
    HZnew = sparse(double(HZ));
    [nX,n0] = size(HX);
    nCurrent = n0;

    rowSupp = cell(0,1);
    rowVals = cell(0,1);
    sourceRow = zeros(0,1);
    gaugeBlocks = struct('sourceXRow',{},'oldSupport',{},'newQubits',{});

    for r = 1:nX
        supp = find(HX(r,:) ~= 0);
        vals = full(HX(r,supp));
        w = numel(supp);

        if w <= targetWeight
            rowSupp{end+1,1} = supp; %#ok<AGROW>
            rowVals{end+1,1} = vals; %#ok<AGROW>
            sourceRow(end+1,1) = r; %#ok<AGROW>
            continue;
        end

        g = w - 3;
        oldN = nCurrent;
        newCols = oldN + (1:g);
        nCurrent = oldN + g;
        HZnew(:,nCurrent) = 0;  % resize sparse matrix

        Zsub = full(HZnew(:,supp));
        Y = zeros(size(HZnew,1),g);

        % First gadget row: x1*z1 + x2*z2 + y1 = 0.
        Y(:,1) = -(vals(1)*Zsub(:,1) + vals(2)*Zsub(:,2));

        % Middle rows: x_{k+2}*z - y_k + y_{k+1} = 0.
        for k = 2:g
            Y(:,k) = Y(:,k-1) - vals(k+1)*Zsub(:,k+1);
        end

        HZnew(:,newCols) = sparse(Y);

        % First row.
        rowSupp{end+1,1} = [supp(1), supp(2), newCols(1)]; %#ok<AGROW>
        rowVals{end+1,1} = [vals(1), vals(2), +1]; %#ok<AGROW>
        sourceRow(end+1,1) = r; %#ok<AGROW>

        % Middle rows.
        for k = 2:g
            rowSupp{end+1,1} = [supp(k+1), newCols(k-1), newCols(k)]; %#ok<AGROW>
            rowVals{end+1,1} = [vals(k+1), -1, +1]; %#ok<AGROW>
            sourceRow(end+1,1) = r; %#ok<AGROW>
        end

        % Last row.
        rowSupp{end+1,1} = [supp(end-1), supp(end), newCols(end)]; %#ok<AGROW>
        rowVals{end+1,1} = [vals(end-1), vals(end), -1]; %#ok<AGROW>
        sourceRow(end+1,1) = r; %#ok<AGROW>

        block.sourceXRow = r;
        block.oldSupport = supp;
        block.newQubits = newCols;
        gaugeBlocks(end+1) = block; %#ok<AGROW>
    end

    % Assemble sparse HX.
    nnzEstimate = sum(cellfun(@numel,rowSupp));
    rr = zeros(nnzEstimate,1);
    cc = zeros(nnzEstimate,1);
    vv = zeros(nnzEstimate,1);
    ptr = 0;
    for r = 1:numel(rowSupp)
        len = numel(rowSupp{r});
        idx = ptr + (1:len);
        rr(idx) = r;
        cc(idx) = rowSupp{r};
        vv(idx) = rowVals{r};
        ptr = ptr + len;
    end
    HXnew = sparse(rr,cc,vv,numel(rowSupp),nCurrent);

    meta.targetWeight = targetWeight;
    meta.numGaugeQubits = nCurrent - n0;
    meta.sourceXRowOfNewRows = sourceRow;
    meta.blocks = gaugeBlocks;
end

%% ========================================================================
%  THICKENING + CHOOSING HEIGHTS OVER Z
%  ========================================================================
function [HXnew, HZnew, meta] = thickenChooseHeightsOverZ(HX, HZ, ell, fixedHeights)
% Signed tensor-product construction.
%
% H is the oriented repetition parity-check matrix with rows [-1,+1].
% The nonbinary tensor-product sign appears in the lower-left Z block:
%
% HXnew = [ HX kron I_ell,  I_nX kron H.' ]
%
% HZfull = [ HZ kron I_ell,                   0
%            -I_n kron H,       HX.' kron I_(ell-1) ]
%
% Then one of the ell copies of every row in the top block is retained.

    HX = sparse(double(HX));
    HZ = sparse(double(HZ));
    [nX,n] = size(HX);
    nZ = size(HZ,1);

    H = orientedRepetitionH(ell);
    Iell = speye(ell);
    Iellm1 = speye(ell-1);

    HXleft = kron(HX, Iell);
    HXright = kron(speye(nX), H.');
    HXnew = sparse([HXleft, HXright]);

    top = sparse([kron(HZ,Iell), sparse(nZ*ell,nX*(ell-1))]);
    bottom = sparse([-kron(speye(n),H), kron(HX.',Iellm1)]);

    % Baseline column load from the bottom block on old-qubit-height columns.
    baseline = full(sum(bottom(:,1:n*ell) ~= 0,1));
    if nargin >= 4 && ~isempty(fixedHeights)
        heights = double(fixedHeights(:));
        if numel(heights) ~= nZ || any(heights < 1) || any(heights > ell) ...
                || any(heights ~= round(heights))
            error('fixedHeights must contain one integer in 1:ell per Z row.');
        end
        loadAfter = reshape(baseline, ell, n).';
        heightCounts = accumarray(heights,1,[ell,1]);
        for j = 1:nZ
            supp = find(HZ(j,:) ~= 0);
            loadAfter(supp,heights(j)) = loadAfter(supp,heights(j)) + 1;
        end
    else
        [heights, loadAfter, heightCounts] = greedyChooseHeights(HZ, ell, baseline);
    end

    keepTop = zeros(nZ,1);
    for j = 1:nZ
        keepTop(j) = (j-1)*ell + heights(j);
    end

    topSelected = top(keepTop,:);
    HZnew = sparse([topSelected; bottom]);

    meta.ell = ell;
    meta.H = H;
    meta.heights = heights;
    meta.heightCounts = heightCounts;
    meta.h = max(heightCounts);
    meta.loadAfter = loadAfter;
    meta.selectedTopOriginalRows = (1:nZ).';
    meta.selectedTopRowIndices = (1:nZ).'; % first rows of HZnew
    meta.bottomRowIndices = (nZ+1:size(HZnew,1)).';
end

function H = orientedRepetitionH(ell)
    rows = zeros(2*(ell-1),1);
    cols = zeros(2*(ell-1),1);
    vals = zeros(2*(ell-1),1);
    ptr = 0;
    for r = 1:ell-1
        rows(ptr+1:ptr+2) = r;
        cols(ptr+1:ptr+2) = [r, r+1];
        vals(ptr+1:ptr+2) = [-1, +1];
        ptr = ptr + 2;
    end
    H = sparse(rows,cols,vals,ell-1,ell);
end

function [heights, load, count] = greedyChooseHeights(HZ, ell, baseline)
% Deterministic greedy choice.  Primary objective: minimize maximum qZ.
% Secondary objectives: balance global height usage and local peak load.

    HZ = sparse(HZ);
    [nZ,n] = size(HZ);
    load = reshape(baseline, ell, n).'; % n x ell
    heights = zeros(nZ,1);
    count = zeros(ell,1);

    rowWeights = full(sum(HZ ~= 0,2));
    [~,order] = sortrows([-rowWeights, (1:nZ).'], [1 2]);

    for idx = 1:nZ
        j = order(idx);
        supp = find(HZ(j,:) ~= 0);
        bestScore = [];
        bestHeight = 1;

        for h = 1:ell
            localNew = load(supp,h) + 1;
            candidateMax = max([load(:); localNew]);
            localPeak = max(localNew);
            localSquare = sum(localNew.^2);
            score = [candidateMax, count(h), localPeak, localSquare, h];

            if isempty(bestScore) || lexicographicLess(score,bestScore)
                bestScore = score;
                bestHeight = h;
            end
        end

        heights(j) = bestHeight;
        load(supp,bestHeight) = load(supp,bestHeight) + 1;
        count(bestHeight) = count(bestHeight) + 1;
    end
end

function tf = lexicographicLess(a,b)
    tf = false;
    for k = 1:numel(a)
        if a(k) < b(k)
            tf = true;
            return;
        elseif a(k) > b(k)
            return;
        end
    end
end

function [tf,ell,heights] = certifiedConditionalBB72Heights(HX,HZ,targetQZ)
% Certified plan after skipping copying and applying gauging to BB72.
    tf = size(HX,1)==144 && size(HX,2)==180 && ...
        size(HZ,1)==36 && size(HZ,2)==180 && targetQZ==3;
    ell = 7;
    heights = [7 7 3 7 7 3 1 6 1 1 6 1 5 4 2 5 4 2 ...
               7 7 3 7 7 3 1 6 1 1 6 1 5 4 2 5 4 2].';
end

function [ell,heights,info] = findMinimalHeightPlanOverZ(HX,HZ,targetQZ,maxEll,maxTime)
% Exact binary assignment MILP for choosing one height per Z row.
% This optimizes qZ; it does not pretend to know dX.
    if max(full(sum(HX~=0,2))) > targetQZ
        error(['No choice of heights can reach target qZ=%d because the ' ...
               'right thickening block already has column weight wX=%d.'], ...
               targetQZ,max(full(sum(HX~=0,2))));
    end

    info.method = 'intlinprog height-feasibility search';
    info.testedEll = [];
    info.exitflags = [];

    if exist('intlinprog','file') ~= 2
        error(['intlinprog is unavailable and no certified plan matches this ' ...
               'input. Supply fixed heights or install Optimization Toolbox.']);
    end

    for ellTry = 2:maxEll
        [h,exitflag] = solveHeightFeasibilityMILP(HZ,ellTry,targetQZ,maxTime);
        info.testedEll(end+1) = ellTry;
        info.exitflags(end+1) = exitflag;
        if ~isempty(h)
            ell = ellTry;
            heights = h;
            info.feasible = true;
            return;
        end
    end
    error('No height assignment found for qZ<=%d with ell<=%d.',targetQZ,maxEll);
end

function [heights,exitflag] = solveHeightFeasibilityMILP(HZ,ell,targetQZ,maxTime)
    HZ = sparse(HZ~=0);
    [nZ,n] = size(HZ);
    Nvar = nZ*ell;

    % Each Z row chooses exactly one height.
    Aeq = sparse(nZ,Nvar);
    for j=1:nZ
        Aeq(j,(j-1)*ell+(1:ell)) = 1;
    end
    beq = ones(nZ,1);

    % In -I kron H, endpoint heights have baseline 1 and internal
    % heights baseline 2.  Top-row load must fit below targetQZ.
    cap = targetQZ-2*ones(n,ell);
    cap(:,1)=targetQZ-1;
    cap(:,ell)=targetQZ-1;

    rowCount = nnz(HZ)*ell;
    rr=zeros(rowCount,1); cc=zeros(rowCount,1); vv=ones(rowCount,1);
    b=zeros(n*ell,1); ptr=0;
    for q=1:n
        zRows=find(HZ(:,q));
        for h=1:ell
            row=(q-1)*ell+h;
            b(row)=cap(q,h);
            len=numel(zRows);
            ids=ptr+(1:len);
            rr(ids)=row;
            cc(ids)=(zRows-1)*ell+h;
            ptr=ptr+len;
        end
    end
    rr=rr(1:ptr); cc=cc(1:ptr); vv=vv(1:ptr);
    A=sparse(rr,cc,vv,n*ell,Nvar);

    f=zeros(Nvar,1);
    intcon=1:Nvar;
    lb=zeros(Nvar,1);
    ub=ones(Nvar,1);
    options=optimoptions('intlinprog','Display','off','MaxTime',maxTime);
    [x,~,exitflag]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
    if isempty(x) || exitflag<=0
        heights=[];
    else
        X=reshape(round(x),ell,nZ).';
        [~,heights]=max(X,[],2);
    end
end

%% ========================================================================
%  SIMULTANEOUS CONING + CELLULATION OVER Z
%  ========================================================================
function [HXnew, HZnew, allData] = coneManyZOverZCellulated( ...
    HX,HZ,zRows,maxCycleLen,targetVertexDegree,maxTries,verbose)
% zRows refer to rows of the input HZ.  All selected rows are removed and
% replaced simultaneously using the block mapping cone.

    HX = sparse(double(HX));
    HZ = sparse(double(HZ));
    zRows = unique(zRows(:).','stable');
    assertCommutesZ(HX,HZ,'input to coning');

    nX = size(HX,1);
    n = size(HX,2);
    nZ = size(HZ,1);

    if any(zRows < 1) || any(zRows > nZ)
        error('Invalid row in zRows.');
    end

    D0list = cell(1,numel(zRows));
    F0list = cell(1,numel(zRows));
    D1list = cell(1,numel(zRows));
    F1list = cell(1,numel(zRows));
    edgeSizes = zeros(1,numel(zRows));
    cycleSizes = zeros(1,numel(zRows));
    supportSizes = zeros(1,numel(zRows));
    allData = repmat(struct(),1,numel(zRows));

    for t = 1:numel(zRows)
        zi = zRows(t);
        [support,zSigns] = coneSupportZ(HZ,zi);
        F1 = coneF1Signed(HZ,zi);
        edges0 = coneEdgesSignedOverZ(HX,HZ,zi);
        D1pre = coneD1Signed(HZ,zi,edges0);
        F0pre = coneF0FromEdges(nX,edges0);
        D0pre = coneD0FundamentalSigned(D1pre);

        checkConeCore(HX,F1,F0pre,D1pre);
        checkConeCycles(D0pre,D1pre);

        [D1,D0,F0,edges,cellInfo] = coneCellulateCyclesBalancedZ( ...
            D1pre,F0pre,edges0,maxCycleLen,targetVertexDegree, ...
            maxTries,1000+zi);

        checkConeCore(HX,F1,F0,D1);
        checkConeCycles(D0,D1);

        D0list{t} = D0;
        F0list{t} = F0;
        D1list{t} = D1;
        F1list{t} = F1;
        edgeSizes(t) = size(D1,1);
        cycleSizes(t) = size(D0,1);
        supportSizes(t) = numel(support);

        allData(t).inputZRow = zi;
        allData(t).support = support;
        allData(t).zSigns = zSigns;
        allData(t).F1 = F1;
        allData(t).F0 = F0;
        allData(t).D1 = D1;
        allData(t).D0 = D0;
        allData(t).F0_pre = F0pre;
        allData(t).D1_pre = D1pre;
        allData(t).D0_pre = D0pre;
        allData(t).edges = edges;
        allData(t).edges_pre = edges0;
        allData(t).cellInfo = cellInfo;
        allData(t).numSupportVertices = numel(support);
        allData(t).numEdges = size(D1,1);
        allData(t).numCycles = size(D0,1);
        allData(t).numCellulationEdges = cellInfo.numAddedEdges;

        if verbose
            fprintf(['  coning row %d/%d (HZ row %d): support=%d, ' ...
                     'edges=%d, cycles=%d, cellulation edges=%d\n'], ...
                t,numel(zRows),zi,numel(support),size(D1,1),size(D0,1), ...
                cellInfo.numAddedEdges);
        end
    end

    Etotal = sum(edgeSizes);
    Ctotal = sum(cycleSizes);
    Wtotal = sum(supportSizes);

    D0block = sparse(Ctotal,Etotal);
    F0all = sparse(nX,Etotal);
    HZtop = sparse(Wtotal,Etotal+n);

    eOff = 0;
    cOff = 0;
    zOff = 0;
    for t = 1:numel(zRows)
        D0 = D0list{t};
        F0 = F0list{t};
        D1 = D1list{t};
        F1 = F1list{t};
        E = size(D1,1);
        C = size(D0,1);
        W = size(D1,2);

        if C > 0
            D0block(cOff+(1:C),eOff+(1:E)) = D0;
        end
        F0all(:,eOff+(1:E)) = F0;
        HZtop(zOff+(1:W),eOff+(1:E)) = D1.';
        HZtop(zOff+(1:W),Etotal+(1:n)) = F1.';

        eOff = eOff + E;
        cOff = cOff + C;
        zOff = zOff + W;
    end

    keep = true(nZ,1);
    keep(zRows) = false;
    HZremaining = HZ(keep,:);

    HXnew = sparse([D0block, sparse(Ctotal,n); F0all, HX]);
    HZnew = sparse([HZtop; sparse(size(HZremaining,1),Etotal), HZremaining]);
    assertCommutesZ(HXnew,HZnew,'simultaneous coning output');
end

function [support,zSigns] = coneSupportZ(HZ,zi)
    if zi < 1 || zi > size(HZ,1)
        error('zi is out of range.');
    end
    support = find(HZ(zi,:) ~= 0);
    zSigns = full(HZ(zi,support));
    if isempty(support)
        error('Cannot cone an empty Z row.');
    end
end

function F1 = coneF1Signed(HZ,zi)
% f1 carries the original signed HZ entries.
    [support,zSigns] = coneSupportZ(HZ,zi);
    F1 = sparse(support,1:numel(support),zSigns,size(HZ,2),numel(support));
end

function edges = coneEdgesSignedOverZ(HX,HZ,zi)
% For each old X row, decompose the signed overlap into flows from positive
% to negative contributions.  Each flow becomes one Xi edge.

    [support,~] = coneSupportZ(HZ,zi);
    z = full(HZ(zi,support));

    owner = zeros(0,1);
    uLocal = zeros(0,1);
    vLocal = zeros(0,1);
    amount = zeros(0,1);

    for r = 1:size(HX,1)
        b = full(HX(r,support)).*z;
        if all(b == 0)
            continue;
        end
        if sum(b) ~= 0
            error('X row %d has nonzero integer overlap with Z row %d.',r,zi);
        end

        pos = find(b > 0);
        neg = find(b < 0);
        posAmt = b(pos);
        negAmt = -b(neg);
        ip = 1;
        ineg = 1;

        while ip <= numel(pos) && ineg <= numel(neg)
            a = min(posAmt(ip),negAmt(ineg));
            owner(end+1,1) = r; %#ok<AGROW>
            uLocal(end+1,1) = pos(ip); %#ok<AGROW>
            vLocal(end+1,1) = neg(ineg); %#ok<AGROW>
            amount(end+1,1) = a; %#ok<AGROW>
            posAmt(ip) = posAmt(ip)-a;
            negAmt(ineg) = negAmt(ineg)-a;
            if posAmt(ip)==0, ip=ip+1; end
            if negAmt(ineg)==0, ineg=ineg+1; end
        end
    end

    edges.ownerXRow = owner;
    edges.uLocal = uLocal;
    edges.vLocal = vLocal;
    edges.uGlobal = support(uLocal).';
    edges.vGlobal = support(vLocal).';
    edges.amount = amount;
    edges.isCellulation = false(numel(owner),1);
    edges.numEdges = numel(owner);
    edges.support = support;
end

function D1 = coneD1Signed(HZ,zi,edges)
% D1 already contains the mapping-cone minus.
% Each edge from positive-overlap vertex u to negative-overlap vertex v has
% D1(e,u)=-a and D1(e,v)=+a.

    [support,~] = coneSupportZ(HZ,zi);
    E = edges.numEdges;
    rr = repelem((1:E).',2);
    cc = reshape([edges.uLocal.'; edges.vLocal.'],[],1);
    vv = reshape([-edges.amount.'; edges.amount.'],[],1);
    D1 = sparse(rr,cc,vv,E,numel(support));
end

function F0 = coneF0FromEdges(nX,edges)
    E = edges.numEdges;
    valid = edges.ownerXRow > 0;
    F0 = sparse(edges.ownerXRow(valid),find(valid),1,nX,E);
end

function D0 = coneD0FundamentalSigned(D1)
    cycles = coneFundamentalCycles(D1);
    D0dense = zeros(numel(cycles),size(D1,1));
    for c = 1:numel(cycles)
        D0dense(c,:) = signedCycleRow(D1,cycles(c).vertices,cycles(c).edges);
    end
    D0 = sparse(D0dense);
end

function [D1cell,D0cell,F0cell,edgesCell,info] = ...
    coneCellulateCyclesBalancedZ(D1pre,F0pre,edgesPre,maxCycleLen, ...
    targetVertexDegree,maxTries,seed)
% Balanced ear cellulation over Z.  It searches over cycle order and ear
% choices so that every D0 row has weight <=4 and every graph vertex has
% final degree <=targetVertexDegree.  Replacement Z rows then have weight
% at most targetVertexDegree+1 because F1 contributes one old-qubit entry.

    if maxCycleLen ~= 4
        error('Balanced cellulation currently requires maxCycleLen=4.');
    end

    cycles = coneFundamentalCycles(D1pre);
    V = size(D1pre,2);
    baseDegree = full(sum(D1pre~=0,1)).';
    capacity = targetVertexDegree-baseDegree;
    if any(capacity<0)
        error('Input coning graph already exceeds target vertex degree.');
    end

    oldState=rng;
    cleanup=onCleanup(@()rng(oldState)); %#ok<NASGU>
    rng(seed,'twister');

    found=false;
    bestPlan=[];
    for attempt=1:maxTries
        order=randperm(numel(cycles));
        load=zeros(V,1);
        plan=cell(numel(cycles),1);
        ok=true;

        for k=1:numel(order)
            c=order(k);
            [faces,chords,loadNew,success] = planCycleEarsBalanced( ...
                cycles(c).vertices,load,capacity);
            if ~success
                ok=false;
                break;
            end
            plan{c}.faces=faces;
            plan{c}.chords=chords;
            load=loadNew;
        end

        if ok && max(baseDegree+load)<=targetVertexDegree
            found=true;
            bestPlan=plan;
            break;
        end
    end

    if ~found
        error(['Balanced cellulation failed after %d attempts. Increase ' ...
               'maxCellulationTries or targetCellulatedVertexDegree.'],maxTries);
    end

    D1cell=full(D1pre);
    F0cell=full(F0pre);
    edgesCell=edgesPre;
    D0dense=zeros(0,size(D1cell,1));
    info.numAddedEdges=0;
    info.addedEdges=struct('edgeIndex',{},'tailLocal',{},'headLocal',{});
    info.attempts=attempt;
    info.targetVertexDegree=targetVertexDegree;

    for c=1:numel(cycles)
        verts=cycles(c).vertices(:).';
        eids=cycles(c).edges(:).';
        L=numel(verts);
        edgeMap=containers.Map('KeyType','char','ValueType','double');

        for k=1:L
            a=verts(k);
            b=verts(mod(k,L)+1);
            edgeMap(edgeKey(a,b))=eids(k);
        end

        chords=bestPlan{c}.chords;
        for k=1:size(chords,1)
            a=chords(k,1);
            b=chords(k,2);
            [D1cell,F0cell,edgesCell,newEdge] = addCellulationEdge( ...
                D1cell,F0cell,edgesCell,a,b);
            D0dense(:,end+1)=0;
            edgeMap(edgeKey(a,b))=newEdge;

            info.numAddedEdges=info.numAddedEdges+1;
            added.edgeIndex=newEdge;
            added.tailLocal=a;
            added.headLocal=b;
            info.addedEdges(end+1)=added; %#ok<AGROW>
        end

        faces=bestPlan{c}.faces;
        for f=1:numel(faces)
            face=faces{f}(:).';
            faceEdges=zeros(1,numel(face));
            for k=1:numel(face)
                a=face(k);
                b=face(mod(k,numel(face))+1);
                key=edgeKey(a,b);
                if ~isKey(edgeMap,key)
                    error('Missing boundary/chord edge while assembling cellulation.');
                end
                faceEdges(k)=edgeMap(key);
            end
            D0dense(end+1,:)=signedCycleRow(sparse(D1cell),face,faceEdges); %#ok<AGROW>
        end
    end

    D1cell=sparse(D1cell);
    D0cell=sparse(D0dense);
    F0cell=sparse(F0cell);
    checkConeCycles(D0cell,D1cell);

    degreeAfter=full(sum(D1cell~=0,1)).';
    info.maxVertexDegree=max(degreeAfter);
    info.degreeAfter=degreeAfter;
    if info.maxVertexDegree>targetVertexDegree
        error('Internal error: cellulated vertex degree exceeded target.');
    end
    if max(full(sum(D0cell~=0,2)))>maxCycleLen
        error('Internal error: a cellulation cycle exceeds maxCycleLen.');
    end
end

function [faces,chords,load,success] = planCycleEarsBalanced(vertices,load,capacity)
% Minimal-chord ear decomposition into triangles/quadrilaterals.
    rem=vertices(:).';
    faces=cell(0,1);
    chords=zeros(0,2);
    success=true;

    while numel(rem)>4
        L=numel(rem);
        candidates=struct('score',{},'face',{},'chord',{},'start',{},'kind',{});

        for i=1:L
            if L==5
                ids=mod((i-1)+(0:2),L)+1;
                face=rem(ids);
                kind=3;
            else
                ids=mod((i-1)+(0:3),L)+1;
                face=rem(ids);
                kind=4;
            end
            a=face(1);
            b=face(end);
            if load(a)+1<=capacity(a) && load(b)+1<=capacity(b)
                candidate.score=[max(load(a)+1,load(b)+1), ...
                    load(a)+load(b),rand];
                candidate.face=face;
                candidate.chord=[a,b];
                candidate.start=i;
                candidate.kind=kind;
                candidates(end+1)=candidate; %#ok<AGROW>
            end
        end

        if isempty(candidates)
            success=false;
            return;
        end

        scores=vertcat(candidates.score);
        [~,ord]=sortrows(scores,[1 2 3]);
        top=min(5,numel(ord));
        chosen=ord(randi(top));
        c=candidates(chosen);
        faces{end+1,1}=c.face; %#ok<AGROW>
        chords(end+1,:)=c.chord; %#ok<AGROW>
        load(c.chord)=load(c.chord)+1;

        rot=rem(mod((c.start-1)+(0:L-1),L)+1);
        if c.kind==4
            rem=[rot(1),rot(4:end)];
        else
            rem=[rot(1),rot(3:end)];
        end
    end

    faces{end+1,1}=rem;
end

function key=edgeKey(a,b)
    a0=min(a,b);
    b0=max(a,b);
    key=sprintf('%d_%d',a0,b0);
end

function [D1cell,D0cell,F0cell,edgesCell,info] = coneCellulateCyclesZ(D1pre,F0pre,edgesPre,maxCycleLen)
% Simple deterministic fan-style cellulation.  It is correct over Z and
% guarantees D0*D1=0, but it is not claimed to minimize wX or wZ.

    D1cell = full(D1pre);
    F0cell = full(F0pre);
    edgesCell = edgesPre;
    cycles = coneFundamentalCycles(D1pre);
    D0dense = zeros(0,size(D1cell,1));
    info.numAddedEdges = 0;
    info.addedEdges = struct('edgeIndex',{},'tailLocal',{},'headLocal',{});

    for c = 1:numel(cycles)
        verts = cycles(c).vertices(:).';
        eids = cycles(c).edges(:).';

        while numel(eids) > maxCycleLen
            split = maxCycleLen;
            v1 = verts(1);
            vk = verts(split);

            [D1cell,F0cell,edgesCell,newEdge] = addCellulationEdge( ...
                D1cell,F0cell,edgesCell,vk,v1);
            D0dense(:,end+1) = 0;

            smallVerts = verts(1:split);
            smallEdges = [eids(1:split-1),newEdge];
            D0dense(end+1,:) = signedCycleRow(sparse(D1cell),smallVerts,smallEdges); %#ok<AGROW>

            verts = [verts(1),verts(split:end)];
            eids = [newEdge,eids(split:end)];

            info.numAddedEdges = info.numAddedEdges + 1;
            added.edgeIndex = newEdge;
            added.tailLocal = vk;
            added.headLocal = v1;
            info.addedEdges(end+1) = added; %#ok<AGROW>
        end

        if ~isempty(eids)
            D0dense(end+1,:) = signedCycleRow(sparse(D1cell),verts,eids); %#ok<AGROW>
        end
    end

    D1cell = sparse(D1cell);
    D0cell = sparse(D0dense);
    F0cell = sparse(F0cell);
    checkConeCycles(D0cell,D1cell);
end

function [D1new,F0new,edgesNew,newEdge] = addCellulationEdge(D1old,F0old,edgesOld,tail,head)
    D1new = D1old;
    newEdge = size(D1old,1)+1;
    newRow = zeros(1,size(D1old,2));
    newRow(tail) = -1;
    newRow(head) = +1;
    D1new(newEdge,:) = newRow;

    F0new = F0old;
    F0new(:,newEdge) = 0;

    edgesNew = edgesOld;
    edgesNew.ownerXRow(newEdge,1) = 0;
    edgesNew.uLocal(newEdge,1) = tail;
    edgesNew.vLocal(newEdge,1) = head;
    edgesNew.uGlobal(newEdge,1) = edgesNew.support(tail);
    edgesNew.vGlobal(newEdge,1) = edgesNew.support(head);
    edgesNew.amount(newEdge,1) = 1;
    edgesNew.isCellulation(newEdge,1) = true;
    edgesNew.numEdges = newEdge;
end

function cycles = coneFundamentalCycles(D1)
% Fundamental cycle basis of the undirected graph represented by D1.

    D1 = sparse(D1);
    E = size(D1,1);
    V = size(D1,2);
    parent = 1:V;
    tree = zeros(0,1);
    nonTree = zeros(0,1);

    for e = 1:E
        [u,v] = edgeEndpoints(D1,e);
        ru = ufFind(parent,u);
        rv = ufFind(parent,v);
        if ru ~= rv
            parent = ufUnion(parent,ru,rv);
            tree(end+1,1) = e; %#ok<AGROW>
        else
            nonTree(end+1,1) = e; %#ok<AGROW>
        end
    end

    adj = cell(1,V);
    for k = 1:numel(tree)
        e = tree(k);
        [u,v] = edgeEndpoints(D1,e);
        adj{u} = [adj{u};v,e]; %#ok<AGROW>
        adj{v} = [adj{v};u,e]; %#ok<AGROW>
    end

    cycles = struct('vertices',{},'edges',{});
    for k = 1:numel(nonTree)
        e = nonTree(k);
        [u,v] = edgeEndpoints(D1,e);
        [pathV,pathE] = treePath(adj,u,v);
        cycles(end+1).vertices = pathV; %#ok<AGROW>
        cycles(end).edges = [pathE,e];
    end
end

function [pathV,pathE] = treePath(adj,startV,targetV)
    V = numel(adj);
    visited = false(1,V);
    prevV = zeros(1,V);
    prevE = zeros(1,V);
    queue = startV;
    visited(startV) = true;

    while ~isempty(queue)
        u = queue(1);
        queue(1) = [];
        if u == targetV, break; end
        nbrs = adj{u};
        for k = 1:size(nbrs,1)
            v = nbrs(k,1);
            e = nbrs(k,2);
            if ~visited(v)
                visited(v) = true;
                prevV(v) = u;
                prevE(v) = e;
                queue(end+1) = v; %#ok<AGROW>
            end
        end
    end

    if ~visited(targetV)
        error('No tree path found.');
    end

    vertsRev = targetV;
    edgesRev = zeros(1,0);
    cur = targetV;
    while cur ~= startV
        edgesRev(end+1) = prevE(cur); %#ok<AGROW>
        cur = prevV(cur);
        vertsRev(end+1) = cur; %#ok<AGROW>
    end
    pathV = fliplr(vertsRev);
    pathE = fliplr(edgesRev);
end

function row = signedCycleRow(D1,vertices,edgeIDs)
    row = zeros(1,size(D1,1));
    L = numel(edgeIDs);
    if numel(vertices) ~= L
        error('Closed cycle representation must have equal vertex/edge counts.');
    end

    for k = 1:L
        a = vertices(k);
        if k < L, b = vertices(k+1); else, b = vertices(1); end
        e = edgeIDs(k);
        [tail,head,amount] = edgeOrientation(D1,e);
        if amount ~= 1
            error(['Cellulation currently requires unit signed incidence. ' ...
                   'Edge %d has magnitude %g.'],e,amount);
        end
        if tail==a && head==b
            row(e)=+1;
        elseif tail==b && head==a
            row(e)=-1;
        else
            error('Edge %d does not connect consecutive cycle vertices.',e);
        end
    end

    if nnz(sparse(row)*D1) ~= 0
        error('Constructed cycle row does not satisfy D0*D1=0.');
    end
end

function [u,v] = edgeEndpoints(D1,e)
    idx = find(D1(e,:) ~= 0);
    if numel(idx) ~= 2
        error('D1 row %d must have exactly two endpoints.',e);
    end
    u = idx(1); v = idx(2);
end

function [tail,head,amount] = edgeOrientation(D1,e)
    idx = find(D1(e,:) ~= 0);
    vals = full(D1(e,idx));
    neg = find(vals<0);
    pos = find(vals>0);
    if numel(neg)~=1 || numel(pos)~=1
        error('D1 row %d must have one negative and one positive endpoint.',e);
    end
    if abs(vals(neg)) ~= abs(vals(pos))
        error('D1 row %d has unequal endpoint magnitudes.',e);
    end
    tail = idx(neg);
    head = idx(pos);
    amount = abs(vals(neg));
end

function root = ufFind(parent,x)
    while parent(x) ~= x
        x = parent(x);
    end
    root = x;
end

function parent = ufUnion(parent,a,b)
    ra = ufFind(parent,a);
    rb = ufFind(parent,b);
    if ra ~= rb
        parent(rb) = ra;
    end
end

function checkConeCore(HX,F1,F0,D1)
    C = F0*D1 + HX*F1;
    if nnz(C) ~= 0
        error('Coning core failed: F0*D1 + HX*F1 ~= 0 over Z.');
    end
end

function checkConeCycles(D0,D1)
    C = D0*D1;
    if nnz(C) ~= 0
        error('Coning cycle condition failed: D0*D1 ~= 0 over Z.');
    end
end

%% ========================================================================
%  STATISTICS, RANK, AND DISTANCE
%  ========================================================================
function s = codeStatsZ(HX,HZ,label)
    s.label = label;
    s.n = size(HX,2);
    s.nX = size(HX,1);
    s.nZ = size(HZ,1);
    s.wX = maxOrZero(full(sum(HX ~= 0,2)));
    s.wZ = maxOrZero(full(sum(HZ ~= 0,2)));
    s.qX = maxOrZero(full(sum(HX ~= 0,1)));
    s.qZ = maxOrZero(full(sum(HZ ~= 0,1)));
    s.minHX = minOrZero(nonzeros(HX));
    s.maxHX = maxOrZero(nonzeros(HX));
    s.minHZ = minOrZero(nonzeros(HZ));
    s.maxHZ = maxOrZero(nonzeros(HZ));
    C = HX*HZ.';
    s.commutatorNNZ = nnz(C);
    s.commutatorMaxAbs = maxAbsSparse(C);
end

function printStats(s,verbose)
    if ~verbose, return; end
    fprintf(['%-27s n=%5d, nX=%5d, nZ=%5d, ' ...
             'wX=%2d, wZ=%2d, qX=%2d, qZ=%2d\n'], ...
        s.label,s.n,s.nX,s.nZ,s.wX,s.wZ,s.qX,s.qZ);
end

function assertCommutesZ(HX,HZ,label)
    if size(HX,2) ~= size(HZ,2)
        error('%s: HX and HZ have different column counts.',label);
    end
    C = HX*HZ.';
    if nnz(C) ~= 0
        error('%s does not commute over Z: nnz=%d, maxAbs=%g.', ...
            label,nnz(C),maxAbsSparse(C));
    end
end

function x = maxAbsSparse(A)
    if nnz(A)==0
        x=0;
    else
        x=max(abs(nonzeros(A)));
    end
end

function x = maxOrZero(v)
    if isempty(v), x=0; else, x=max(v); end
end

function x = minOrZero(v)
    if isempty(v), x=0; else, x=min(v); end
end

function r = rankModPrime(A,p)
% Sparse dictionary-style Gaussian elimination over F_p.
% This is exact but can take time for very large matrices.

    if ~isprime(p)
        error('p must be prime.');
    end
    A = sparse(mod(A,p));

    if p == 2
        r = gf2RankIncremental(A);
        return;
    end

    basis = cell(1,size(A,2));
    hasPivot = false(1,size(A,2));
    r = 0;

    for i = 1:size(A,1)
        row = full(A(i,:));
        row = mod(row,p);
        while any(row)
            pivot = find(row~=0,1,'last');
            if hasPivot(pivot)
                coeff = row(pivot);
                row = mod(row - coeff*basis{pivot},p);
            else
                invPivot = modInverse(row(pivot),p);
                row = mod(invPivot*row,p);
                basis{pivot} = row;
                hasPivot(pivot) = true;
                r = r+1;
                break;
            end
        end
    end
end

function inva = modInverse(a,p)
    a = mod(a,p);
    for x = 1:p-1
        if mod(a*x,p)==1
            inva=x;
            return;
        end
    end
    error('No modular inverse found.');
end

function r = gf2RankIncremental(A)
    A = logical(mod(A,2));
    n = size(A,2);
    basis = cell(1,n);
    hasPivot = false(1,n);
    r = 0;

    for i = 1:size(A,1)
        v = full(A(i,:));
        while any(v)
            p = find(v,1,'first');
            if hasPivot(p)
                v = xor(v,basis{p});
            else
                basis{p}=v;
                hasPivot(p)=true;
                r=r+1;
                break;
            end
        end
    end
end

function D = distanceReport(HX,HZ,opts,copyMeta)
    D = struct();
    D.inputDistance = 6;
    D.field = 'binary/qubit distance uses final matrices reduced mod 2';
    D.theory.copying_dX = 6;
    D.theory.copying_dZ = copyMeta.qXcopies*6;
    D.theory.afterGauging_dZ_lowerBound = D.theory.copying_dZ;
    D.theory.afterGauging_dX = 'not fixed by Lemma 2';
    D.theory.afterThickening = sprintf('dX is multiplied by ell=%d; dZ is unchanged',opts.ell);
    D.theory.afterConing = 'dX does not decrease; no general exact dZ follows from weights';
    D.exact = NaN;
    D.dX = NaN;
    D.dZ = NaN;
    D.status = 'not computed';

    mode = lower(string(opts.distanceMode));
    if mode == "report_only"
        D.status = ['Exact final distance was not run. Set opts.distanceMode=' ...
                    '''milp_binary'' to attempt it with intlinprog.'];
        return;
    elseif mode ~= "milp_binary"
        error('Unknown distanceMode: %s',opts.distanceMode);
    end

    if exist('intlinprog','file') ~= 2
        D.status = ['Optimization Toolbox/intlinprog is unavailable. ' ...
                    'Exact final distance was not computed.'];
        return;
    end

    try
        [dX,dZ,details] = cssDistanceBinaryMILP(HX,HZ,opts.distanceTimePerSolve);
        D.dX = dX;
        D.dZ = dZ;
        D.exact = min(dX,dZ);
        D.details = details;
        if details.allExact
            D.status = 'exact binary distance computed by MILP';
        else
            D.status = ['MILP returned best feasible upper bounds, but at least ' ...
                        'one solve hit its time limit; exactness is not certified.'];
        end
    catch ME
        D.status = ['Distance attempt failed: ' ME.message];
    end
end

function [dX,dZ,details] = cssDistanceBinaryMILP(HX,HZ,timePerSolve)
% Exact when every constituent MILP terminates optimally.
% This can be very expensive for the final 4000+ qubit code.

    HX2 = sparse(mod(HX,2));
    HZ2 = sparse(mod(HZ,2));

    LX = quotientLogicalBasis(HZ2,HX2); % X logicals: ker(HZ)/row(HX)
    LZ = quotientLogicalBasis(HX2,HZ2); % Z logicals: ker(HX)/row(HZ)

    [dZ,zInfo] = minOddLogicalMILP(HX2,LX,timePerSolve);
    [dX,xInfo] = minOddLogicalMILP(HZ2,LZ,timePerSolve);

    details.XlogicalBasis = LX;
    details.ZlogicalBasis = LZ;
    details.Xsolves = xInfo;
    details.Zsolves = zInfo;
    details.allExact = xInfo.allExact && zInfo.allExact;
end

function L = quotientLogicalBasis(Hcomm,Hstab)
% Basis for ker(Hcomm)/row(Hstab) over F2.

    N = gf2Nullspace(Hcomm);
    base = independentRowBasis(Hstab);
    [~,basisRows,pivots] = initializeGF2Basis(base); %#ok<ASGLU>
    logicals = false(0,size(Hcomm,2));

    for i = 1:size(N,1)
        [isNew,basisRows,pivots,reduced] = addToGF2Basis(N(i,:),basisRows,pivots);
        if isNew
            logicals(end+1,:) = reduced; %#ok<AGROW>
        end
    end
    L = sparse(logicals);
end

function B = independentRowBasis(A)
    A = logical(mod(A,2));
    [~,rows,~] = initializeGF2Basis(A);
    if isempty(rows)
        B = sparse(0,size(A,2));
    else
        B = sparse(vertcat(rows{:}));
    end
end

function [rankValue,basisRows,pivots] = initializeGF2Basis(A)
    basisRows = cell(0,1);
    pivots = zeros(0,1);
    for i = 1:size(A,1)
        [~,basisRows,pivots] = addToGF2Basis(full(A(i,:)),basisRows,pivots);
    end
    rankValue = numel(pivots);
end

function [isNew,basisRows,pivots,reduced] = addToGF2Basis(v,basisRows,pivots)
    v = logical(v);
    for k = 1:numel(pivots)
        if v(pivots(k))
            v = xor(v,basisRows{k});
        end
    end
    reduced = v;
    if ~any(v)
        isNew=false;
        return;
    end

    p = find(v,1,'first');
    insertAt = find(pivots>p,1);
    if isempty(insertAt)
        pivots(end+1,1)=p;
        basisRows{end+1,1}=v;
    else
        pivots = [pivots(1:insertAt-1);p;pivots(insertAt:end)];
        basisRows = [basisRows(1:insertAt-1);{v};basisRows(insertAt:end)];
    end
    isNew=true;
end

function N = gf2Nullspace(A)
% Dense RREF-based nullspace.  Used only by optional distance solver.

    A = logical(full(mod(A,2)));
    [m,n] = size(A);
    row = 1;
    pivots = zeros(0,1);

    for col = 1:n
        if row > m, break; end
        idx = find(A(row:end,col),1);
        if isempty(idx), continue; end
        idx = idx+row-1;
        if idx~=row
            temp=A(row,:); A(row,:)=A(idx,:); A(idx,:)=temp;
        end
        other = find(A(:,col));
        other(other==row)=[];
        for q = other.'
            A(q,:) = xor(A(q,:),A(row,:));
        end
        pivots(end+1,1)=col; %#ok<AGROW>
        row=row+1;
    end

    free = setdiff(1:n,pivots);
    N = false(numel(free),n);
    for k = 1:numel(free)
        f = free(k);
        N(k,f)=true;
        for r = 1:numel(pivots)
            if A(r,f)
                N(k,pivots(r))=true;
            end
        end
    end
    N = sparse(N);
end

function [best,info] = minOddLogicalMILP(Hcomm,Lopposite,timePerSolve)
% Minimize weight subject to commuting with Hcomm and having odd overlap
% with at least one basis logical of the opposite type.  We solve one MILP
% per basis logical and take the minimum.

    H = sparse(mod(Hcomm,2));
    L = sparse(mod(Lopposite,2));
    n = size(H,2);
    m = size(H,1);
    k = size(L,1);

    if k==0
        best=Inf;
        info.allExact=true;
        info.objectives=[];
        info.exitflags=[];
        return;
    end

    objectives = Inf(k,1);
    exitflags = zeros(k,1);
    options = optimoptions('intlinprog','Display','off','MaxTime',timePerSolve);

    rowWeights = full(sum(H,2));
    tUpper = floor(rowWeights/2);

    for i = 1:k
        l = L(i,:);
        sUpper = floor(nnz(l)/2);
        Nvar = n+m+1;
        f = [ones(n,1);zeros(m+1,1)];
        intcon = 1:Nvar;
        Aeq = [H,-2*speye(m),sparse(m,1); ...
               l,sparse(1,m),-2];
        beq = [zeros(m,1);1];
        lb = zeros(Nvar,1);
        ub = [ones(n,1);tUpper;sUpper];

        [~,fval,exitflag] = intlinprog(f,intcon,[],[],Aeq,beq,lb,ub,options);
        exitflags(i)=exitflag;
        if ~isempty(fval)
            objectives(i)=round(fval);
        end
    end

    best=min(objectives);
    info.objectives=objectives;
    info.exitflags=exitflags;
    info.allExact=all(exitflags>0);
end
