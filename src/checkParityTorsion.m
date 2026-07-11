function R = checkParityTorsion(HX, HZ)
%CHECKPARITYTORSION Exact torsion test for integer parity-check matrices.
%
%   R = checkParityTorsion(H)
%       Checks torsion in coker(H) = Z^(number of rows)/im(H).
%
%   R = checkParityTorsion(HX,HZ)
%       Checks:
%         1) torsion in coker(HX),
%         2) torsion in coker(HZ),
%         3) torsion in the integral CSS homology groups
%
%              H_Zlogical = ker_Z(HX) / im_Z(HZ.')
%              H_Xlogical = ker_Z(HZ) / im_Z(HX.').
%
%   Homology is computed only if HX*HZ.' = 0 exactly over Z.
%
%   REQUIRES: Symbolic Math Toolbox.
%
%   A Smith invariant factor d_i > 1 contributes a torsion factor
%
%                         Z / d_i Z.
%
%   Example:
%       R = checkParityTorsion(HX, HZ_LDI);
%
%       % Check only one parity-check matrix:
%       R = checkParityTorsion(HX);
%
%   IMPORTANT:
%     This is an exact integer calculation intended for small or moderate
%     matrices. Smith normal form can be very expensive for matrices with
%     thousands of rows/columns. Do not expect the 3278-by-6626 final BB
%     matrices to finish quickly with an exact symbolic calculation.
%
%   The rectangular Smith-normal-form routine is implemented locally;
%   MATLAB's built-in smithForm is not used because it is restricted to
%   square matrices in current documented versions.

    narginchk(1,2);

    validateIntegerMatrix(HX, 'HX');
    HXs = sym(full(HX));

    fprintf('\n============================================================\n');
    fprintf('Exact integer torsion analysis using Smith normal form\n');
    fprintf('============================================================\n');

    R = struct;
    R.cokerHX = smithCokernelReport(HXs, 'coker(HX)');

    if nargin < 2 || isempty(HZ)
        return;
    end

    validateIntegerMatrix(HZ, 'HZ');
    HZs = sym(full(HZ));

    if size(HX,2) ~= size(HZ,2)
        error('HX and HZ must have the same number of columns.');
    end

    R.cokerHZ = smithCokernelReport(HZs, 'coker(HZ)');

    commutator = HXs * HZs.';
    R.integerCommutation = isZeroSym(commutator);

    fprintf('\nIntegral CSS commutation HX*HZ.'' = 0: %s\n', ...
        yesNo(R.integerCommutation));

    if ~R.integerCommutation
        fprintf(['CSS homology torsion was skipped because this pair does not\n' ...
                 'form a chain complex over Z. Use the signed LDI pair for\n' ...
                 'integer homology, or reduce modulo p for field analysis.\n']);
        R.ZlogicalHomology = [];
        R.XlogicalHomology = [];
        return;
    end

    % Z logical classes: ker(HX)/im(HZ.')
    R.ZlogicalHomology = integralHomologyReport( ...
        HXs, HZs.', 'H_Zlogical = ker_Z(HX) / im_Z(HZ.'')');

    % X logical classes: ker(HZ)/im(HX.')
    R.XlogicalHomology = integralHomologyReport( ...
        HZs, HXs.', 'H_Xlogical = ker_Z(HZ) / im_Z(HX.'')');
end


function out = smithCokernelReport(A, label)
% For A : Z^n -> Z^m, report coker(A) = Z^m/im(A).

    m = size(A,1);
    n = size(A,2);

    if m == 0 || n == 0
        out = emptyCokernelReport(m,n);
        printReport(out,label);
        return;
    end

    [D,U,V] = smithNormalFormZ(A);

    if ~isZeroSym(D - U*A*V)
        error('Internal Smith decomposition verification D = U*A*V failed.');
    end

    [nonzeroFactors,torsionFactors,rankA] = extractSmithFactors(D);

    out = struct;
    out.matrixSize = [m,n];
    out.rankOverQ = rankA;
    out.freeRank = m-rankA;
    out.nonzeroInvariantFactors = nonzeroFactors;
    out.torsionInvariantFactors = torsionFactors;
    out.hasTorsion = ~isempty(torsionFactors);
    out.torsionPrimes = extractTorsionPrimes(torsionFactors);

    printReport(out,label);
end


function out = integralHomologyReport(A,B,label)
% Compute ker_Z(A)/im_Z(B), assuming A*B = 0.
% A is m-by-n and B is n-by-s.

    if size(A,2) ~= size(B,1)
        error('For homology, size(A,2) must equal size(B,1).');
    end

    if ~isZeroSym(A*B)
        error('The supplied maps do not satisfy A*B = 0 over Z.');
    end

    n = size(A,2);

    if isempty(A)
        rankA = 0;
        V = sym(eye(n));
    else
        [D,U,V] = smithNormalFormZ(A);
        if ~isZeroSym(D-U*A*V)
            error('Smith verification failed for the outgoing map.');
        end
        [~,~,rankA] = extractSmithFactors(D);
    end

    nullityA = n-rankA;

    % If D = U*A*V, then the last n-rank(A) columns of V form an
    % integral basis for ker_Z(A).  Y = V\B gives coordinates of im(B)
    % in that basis; its first rank(A) rows must be zero.
    if size(B,2) == 0
        C = sym(zeros(nullityA,0));
    else
        Y = simplify(V\B);

        if rankA > 0 && ~isZeroSym(Y(1:rankA,:))
            error('im(B) is not contained in ker(A) in Smith coordinates.');
        end

        C = Y(rankA+1:end,:);
    end

    if isempty(C) || size(C,2) == 0
        nonzeroFactors = sym([]);
        torsionFactors = sym([]);
        rankC = 0;
    else
        [DC,UC,VC] = smithNormalFormZ(C);
        if ~isZeroSym(DC-UC*C*VC)
            error('Smith verification failed for the homology quotient.');
        end
        [nonzeroFactors,torsionFactors,rankC] = extractSmithFactors(DC);
    end

    out = struct;
    out.outgoingMapRankOverQ = rankA;
    out.kernelRank = nullityA;
    out.imageRankInsideKernel = rankC;
    out.freeRank = nullityA-rankC;
    out.nonzeroInvariantFactors = nonzeroFactors;
    out.torsionInvariantFactors = torsionFactors;
    out.hasTorsion = ~isempty(torsionFactors);
    out.torsionPrimes = extractTorsionPrimes(torsionFactors);

    printReport(out,label);
end


function [D,U,V] = smithNormalFormZ(A)
%SMITHNORMALFORMZ Smith normal form of a rectangular integer matrix.
%
% Returns unimodular U,V and diagonal D satisfying D = U*A*V.
% Arithmetic is symbolic/exact. The algorithm repeatedly uses Bezout row
% and column transformations and enforces divisibility of each pivot into
% the remaining submatrix.

    A = sym(A);
    [m,n] = size(A);
    U = sym(eye(m));
    V = sym(eye(n));

    kmax = min(m,n);

    for k = 1:kmax
        % Find any nonzero element in the unreduced submatrix.
        block = A(k:end,k:end);
        mask = ~isAlways(block == 0);
        linearIndex = find(mask,1,'first');

        if isempty(linearIndex)
            break;
        end

        [ii,jj] = ind2sub(size(block),linearIndex);
        pivotRow = k+ii-1;
        pivotCol = k+jj-1;

        if pivotRow ~= k
            A([k,pivotRow],:) = A([pivotRow,k],:);
            U([k,pivotRow],:) = U([pivotRow,k],:);
        end

        if pivotCol ~= k
            A(:,[k,pivotCol]) = A(:,[pivotCol,k]);
            V(:,[k,pivotCol]) = V(:,[pivotCol,k]);
        end

        iteration = 0;
        while true
            iteration = iteration+1;
            if iteration > 100000
                error('Smith algorithm did not converge at pivot %d.',k);
            end

            if isAlways(A(k,k) < 0)
                A(k,:) = -A(k,:);
                U(k,:) = -U(k,:);
            end

            % Reduce/clear entries below the pivot.
            for i = k+1:m
                b = A(i,k);
                if isAlways(b == 0)
                    continue;
                end

                a = A(k,k);
                if isAlways(mod(b,a) == 0)
                    q = b/a;
                    A(i,:) = A(i,:) - q*A(k,:);
                    U(i,:) = U(i,:) - q*U(k,:);
                else
                    [g,s,t] = extendedGCDZ(a,b);
                    M = [s, t; -b/g, a/g];
                    A([k,i],:) = M*A([k,i],:);
                    U([k,i],:) = M*U([k,i],:);
                end
            end

            % Reduce/clear entries to the right of the pivot.
            for j = k+1:n
                b = A(k,j);
                if isAlways(b == 0)
                    continue;
                end

                a = A(k,k);
                if isAlways(mod(b,a) == 0)
                    q = b/a;
                    A(:,j) = A(:,j) - q*A(:,k);
                    V(:,j) = V(:,j) - q*V(:,k);
                else
                    [g,s,t] = extendedGCDZ(a,b);
                    N = [s, -b/g; t, a/g];
                    A(:,[k,j]) = A(:,[k,j])*N;
                    V(:,[k,j]) = V(:,[k,j])*N;
                end
            end

            % Column operations can reintroduce entries below the pivot,
            % and row operations can reintroduce entries to its right.
            if ~isZeroSym(A(k+1:end,k)) || ~isZeroSym(A(k,k+1:end))
                continue;
            end

            pivot = A(k,k);
            offendingRow = 0;
            offendingCol = 0;

            % Smith condition: pivot must divide every entry remaining.
            for i = k+1:m
                for j = k+1:n
                    if ~isAlways(mod(A(i,j),pivot) == 0)
                        offendingRow = i;
                        offendingCol = j;
                        break;
                    end
                end
                if offendingRow ~= 0
                    break;
                end
            end

            if offendingRow == 0
                break;
            end

            % Bring an offending entry into the pivot row. The next column
            % Bezout operation replaces the pivot by gcd(pivot,offender).
            A(k,:) = A(k,:) + A(offendingRow,:);
            U(k,:) = U(k,:) + U(offendingRow,:);

            %#ok<NASGU> offendingCol is retained for conceptual clarity.
        end

        if isAlways(A(k,k) < 0)
            A(k,:) = -A(k,:);
            U(k,:) = -U(k,:);
        end
    end

    D = A;

    % Final structural checks.
    if ~isDiagonalRectangular(D)
        error('The computed Smith form is not diagonal.');
    end

    diagonalLength = min(size(D));
    previous = sym(0);
    for i = 1:diagonalLength
        current = abs(D(i,i));
        if isAlways(current == 0)
            continue;
        end
        if ~isAlways(previous == 0) && ~isAlways(mod(current,previous) == 0)
            error('Smith divisibility condition d_i | d_{i+1} failed.');
        end
        previous = current;
    end
end


function [g,s,t] = extendedGCDZ(a,b)
% Return g > 0 and Bezout coefficients s,t with s*a+t*b = g.

    a = sym(a);
    b = sym(b);

    aa = abs(a);
    bb = abs(b);

    oldR = aa;
    r = bb;
    oldS = sym(1);
    s0 = sym(0);
    oldT = sym(0);
    t0 = sym(1);

    while ~isAlways(r == 0)
        q = floor(oldR/r);

        temp = oldR-q*r;
        oldR = r;
        r = temp;

        temp = oldS-q*s0;
        oldS = s0;
        s0 = temp;

        temp = oldT-q*t0;
        oldT = t0;
        t0 = temp;
    end

    g = oldR;
    s = oldS*sign(a);
    t = oldT*sign(b);

    if ~isAlways(s*a+t*b == g)
        error('Extended-GCD verification failed.');
    end
end


function [nonzeroFactors,torsionFactors,rankA] = extractSmithFactors(D)
    diagonalLength = min(size(D));
    nonzeroFactors = sym([]);
    torsionFactors = sym([]);
    rankA = 0;

    for i = 1:diagonalLength
        di = abs(D(i,i));
        if ~isAlways(di == 0)
            rankA = rankA+1;
            nonzeroFactors(end+1) = di; %#ok<AGROW>
            if ~isAlways(di == 1)
                torsionFactors(end+1) = di; %#ok<AGROW>
            end
        end
    end
end


function primesOut = extractTorsionPrimes(torsionFactors)
    primesOut = [];

    for i = 1:numel(torsionFactors)
        value = double(torsionFactors(i));
        if isfinite(value) && value >= 2 && value <= flintmax
            primesOut = union(primesOut,unique(factor(round(value))));
        end
    end
end


function printReport(out,label)
    fprintf('\n--- %s ---\n',label);

    if isfield(out,'matrixSize')
        fprintf('Matrix size: %d x %d\n',out.matrixSize(1),out.matrixSize(2));
        fprintf('Rank over Q: %d\n',out.rankOverQ);
    else
        fprintf('Kernel rank before quotient: %d\n',out.kernelRank);
        fprintf('Image rank inside kernel: %d\n',out.imageRankInsideKernel);
    end

    fprintf('Free rank: %d\n',out.freeRank);
    fprintf('Torsion present: %s\n',yesNo(out.hasTorsion));

    if out.hasTorsion
        fprintf('Torsion invariant factors d_i > 1:\n');
        disp(out.torsionInvariantFactors(:).');

        if ~isempty(out.torsionPrimes)
            fprintf('Torsion primes detected: ');
            fprintf('%d ',out.torsionPrimes);
            fprintf('\n');
        end

        fprintf('Torsion subgroup = direct sum of Z/d_i Z.\n');
    else
        fprintf('No finite torsion factors were found.\n');
    end
end


function out = emptyCokernelReport(m,n)
    out = struct;
    out.matrixSize = [m,n];
    out.rankOverQ = 0;
    out.freeRank = m;
    out.nonzeroInvariantFactors = sym([]);
    out.torsionInvariantFactors = sym([]);
    out.hasTorsion = false;
    out.torsionPrimes = [];
end


function validateIntegerMatrix(A,name)
    if ~(isnumeric(A) || islogical(A))
        error('%s must be a numeric or logical MATLAB matrix.',name);
    end

    values = full(double(A(:)));

    if any(~isfinite(values))
        error('%s contains NaN or Inf.',name);
    end

    if any(values ~= round(values))
        error('%s must contain integer entries.',name);
    end
end


function tf = isZeroSym(M)
    if isempty(M)
        tf = true;
    else
        tf = all(all(isAlways(M == 0)));
    end
end


function tf = isDiagonalRectangular(M)
    [m,n] = size(M);
    tf = true;

    for i = 1:m
        for j = 1:n
            if i ~= j && ~isAlways(M(i,j) == 0)
                tf = false;
                return;
            end
        end
    end
end


function text = yesNo(tf)
    if tf
        text = 'YES';
    else
        text = 'NO';
    end
end
