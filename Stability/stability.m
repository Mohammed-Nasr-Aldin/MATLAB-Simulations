function routhStability()
%% ===================== Control Systems Suite =====================
%  1) Routh stability   2) Root Locus   3) Bode plot
%  Save this file as routhStability.m, then type  routhStability
clc; close all;

again = true;
while again
    fprintf('\n=========== Control Systems Suite ===========\n');
    fprintf('  1) Routh stability\n');
    fprintf('  2) Root Locus\n');
    fprintf('  3) Bode plot\n');
    choice = readChoice('Choose analysis type (1/2/3): ', [1 2 3]);
    try
        switch choice
            case 1, runRouth();
            case 2, runRootLocus();
            case 3, runBode();
        end
    catch ME
        fprintf('\nSomething went wrong: %s\n', ME.message);
    end
    again = readChoice(sprintf('\nRun another analysis? (1 = yes, 0 = no): '), [0 1]) == 1;
end
fprintf('\nDone.\n');
end

%% ===================================================================
%                              ROUTH
%% ===================================================================
function runRouth()
    fprintf('\n--- Routh-Hurwitz ---\n');
    fprintf('  1) Numeric coefficients\n');
    fprintf('  2) Coefficients containing a gain K\n');
    mode = readChoice('Choose mode (1 or 2): ', [1 2]);
    n = readInt('Enter the order n of the system = ', 1, 30);

    if mode == 1
        fprintf('\nCharacteristic eq: c_n*S^n + ... + c_1*S + c_0 = 0\n');
        while true
            S = zeros(1, n+1);
            for k = 1:n+1
                S(k) = readNum(sprintf('Enter coefficient of S^%d = ', n-k+1));
            end
            if any(S ~= 0), break; end
            disp('All coefficients are zero. Please enter them again.');
        end
        S = trimLead(S);
        n = numel(S) - 1;
        if n < 1
            disp('The polynomial reduces to a constant: there are no poles to analyze.');
            return;
        end
        fprintf('\nActual order of the system = %d\n', n);

        A = analyzePoly(S);
        printRouthTable(A.R, A.n, A.info);
        printAux(A.info, A.n);

        fprintf('\nFirst column: %s\n', mat2str(A.R(:,1)', 4));
        fprintf('\nFrom Routh table:\n');
        fprintf('  Poles in RHP            = %d\n', A.numRHP);
        fprintf('  Poles on imaginary axis = %d\n', A.numImag);
        fprintf('  Poles in LHP            = %d\n', A.numLHP);

        r = roots(S);
        tol = 1e-6;
        repeated = hasRepeatedImag(r);

        if A.numRHP > 0
            resultTxt = 'UNSTABLE';
        elseif A.numImag > 0 && repeated
            resultTxt = 'UNSTABLE (repeated poles on jw axis)';
        elseif A.numImag > 0
            resultTxt = 'MARGINALLY STABLE';
        else
            resultTxt = 'STABLE';
        end
        fprintf('Result: system is %s\n', resultTxt);

        fprintf('\nVerification with roots():\n');
        disp(r);
        nR = sum(real(r) > tol);  nI = sum(abs(real(r)) <= tol);  nL = sum(real(r) < -tol);
        fprintf('  RHP = %d | Imag axis = %d | LHP = %d\n', nR, nI, nL);
        if nR ~= A.numRHP || nI ~= A.numImag
            fprintf('  Warning: roots() disagrees with the Routh counts (numerical sensitivity near the jw axis).\n');
        end
        plotPoles(r, resultTxt, A.numRHP);

    else
        if isempty(ver('symbolic'))
            disp('Mode 2 needs the Symbolic Math Toolbox, which is not installed.');
            return;
        end
        syms K real
        fprintf('\nEnter each coefficient as an expression in K, e.g.  3   or   K+2   or   2*K   (2K also works)\n');
        S = sym(zeros(1, n+1));
        for k = 1:n+1
            while true
                str = strtrim(input(sprintf('Coefficient of S^%d = ', n-k+1), 's'));
                str = strrep(str, 'k', 'K');
                if isempty(str) || isempty(regexp(str, '^[0-9K\.\+\-\*/\^\(\)\s]+$', 'once'))
                    disp('Invalid expression. Use only numbers, K, + - * / ^ and parentheses.');
                    continue;
                end
                str = regexprep(str, '(\d|\))\s*K', '$1*K');      % 2K  -> 2*K
                str = regexprep(str, '(\d|K|\))\s*\(', '$1*(');   % 2(K+1) -> 2*(K+1)
                try
                    S(k) = sym(eval(str));
                    break;
                catch
                    disp('Could not understand that expression. Try again.');
                end
            end
        end

        while numel(S) > 1 && symZero(S(1))
            S(1) = [];
        end
        n = numel(S) - 1;
        if n < 1
            disp('The polynomial reduces to a constant: nothing to analyze.');
            return;
        end
        levelNote = ~isempty(symvar(S(1)));
        if levelNote
            disp('Note: the leading coefficient depends on K; the order drops at the K that makes it zero.');
        end

        % ---------- symbolic Routh table (stops if a pivot is identically zero) ----------
        cols = ceil((n+1)/2) + 1;
        R = sym(zeros(n+1, cols));
        r1 = S(1:2:end);  R(1, 1:numel(r1)) = r1;
        r2 = S(2:2:end);  R(2, 1:numel(r2)) = r2;
        lastRow = n + 1;
        tableOK = true;
        for i = 3:n+1
            piv = R(i-1,1);
            if symZero(piv)
                tableOK = false;  lastRow = i - 1;  break;
            end
            for j = 1:cols-1
                R(i,j) = simplify((piv*R(i-2,j+1) - R(i-2,1)*R(i-1,j+1)) / piv);
            end
        end
        firstCol = simplify(R(1:lastRow, 1));

        fprintf('\nFirst column of the Routh table:\n');
        for i = 1:lastRow
            fprintf('S^%-2d | ', n-i+1);
            disp(firstCol(i));
        end
        if ~tableOK
            fprintf(['\nSpecial case: the pivot of row S^%d is zero for EVERY K (zero row / zero pivot).\n' ...
                     'The symbolic table stops here; stability is checked numerically for each K range,\n' ...
                     'using the epsilon / auxiliary-polynomial rules.\n'], n-lastRow+1);
        end

        % ---------- critical K values ----------
        crit = [];
        % (a) K where a first-column element is zero or its denominator is zero
        for i = 1:lastRow
            try
                [nu, de] = numden(firstCol(i));
                for e = {nu, de}
                    ex = e{1};
                    if ~isempty(symvar(ex))
                        sol = double(solve(ex == 0, K));
                        sol = sol(abs(imag(sol)) < 1e-9);
                        crit = [crit; real(sol(:))]; %#ok<AGROW>
                    end
                end
            catch
            end
        end
        % (b) K where the leading coefficient vanishes
        if levelNote
            try
                sol = double(solve(S(1) == 0, K));
                sol = sol(abs(imag(sol)) < 1e-9);
                crit = [crit; real(sol(:))]; %#ok<AGROW>
            catch
            end
        end
        % (c) K where a root sits on the jw axis (solve P(jw) = 0 for w and K)
        try
            syms w real
            Pjw = expand(sum(S .* (1j*w).^(n:-1:0)));
            sol = solve([real(Pjw) == 0, imag(Pjw) == 0], [w K], 'Real', true);
            if isstruct(sol)
                wv = double(sol.w);  Kv = double(sol.K);
            else
                wv = double(sol(:,1));  Kv = double(sol(:,2));
            end
            okv = isfinite(wv) & isfinite(Kv) & abs(imag(wv)) < 1e-9 & abs(imag(Kv)) < 1e-9;
            crit = [crit; real(Kv(okv))]; %#ok<AGROW>
        catch
        end
        crit = unique(round(crit(isfinite(crit)), 8));

        % ---------- stability by K range (numeric Routh with all special cases) ----------
        edges = [-inf; crit(:); inf];
        fprintf('\nStability by K range:\n');
        stableIdx = [];
        for k = 1:numel(edges)-1
            lo = edges(k); hi = edges(k+1);
            if isinf(lo) && isinf(hi), pt = 0;
            elseif isinf(lo),          pt = hi - 1;
            elseif isinf(hi),          pt = lo + 1;
            else,                      pt = (lo + hi) / 2;
            end
            [isSt, st] = evalAtK(S, K, pt);
            if isSt, stableIdx(end+1) = k; end %#ok<AGROW>
            fprintf('  %s : %s\n', rangeStr(lo, hi), st);
        end

        if isempty(stableIdx)
            fprintf('\nNo value of K makes the system stable.\n');
        else
            fprintf('\nStable for:\n');
            for k = stableIdx
                fprintf('  %s\n', rangeStr(edges(k), edges(k+1)));
            end
        end

        % ---------- behaviour exactly at each critical K ----------
        if ~isempty(crit)
            fprintf('\nCritical K values:\n');
            for kc = crit(:)'
                try
                    c = double(subs(S, K, kc));
                    c = trimLead(c);
                    Ac = analyzePoly(c);
                    rr = roots(c);
                    onAx = rr(abs(real(rr)) <= 1e-6);
                    wj = unique(round(abs(imag(onAx)), 6));
                    if Ac.numRHP > 0
                        fprintf('  K = %.6g -> UNSTABLE (%d RHP poles)\n', kc, Ac.numRHP);
                    elseif Ac.numImag > 0
                        fprintf('  K = %.6g -> MARGINALLY STABLE, poles on jw axis at w = %s rad/s\n', ...
                                kc, mat2str(wj', 5));
                    else
                        fprintf('  K = %.6g -> stable (a first-column element changes sign here)\n', kc);
                    end
                catch
                    fprintf('  K = %.6g (could not be analyzed)\n', kc);
                end
            end
        end
    end
end

function [isStable, txt] = evalAtK(S, K, pt)
    try
        c = double(subs(S, K, pt));
        A = analyzePoly(c);
        if A.numRHP == 0 && A.numImag == 0
            isStable = true;   txt = 'STABLE';
        elseif A.numRHP == 0
            isStable = false;  txt = 'MARGINALLY STABLE';
        else
            isStable = false;  txt = sprintf('UNSTABLE (%d RHP poles)', A.numRHP);
        end
    catch
        isStable = false;  txt = 'could not be evaluated';
    end
end

function tf = symZero(x)
    tf = isequal(simplify(sym(x)), sym(0));
end

%% ===================================================================
%                           ROOT LOCUS
%% ===================================================================
function runRootLocus()
    fprintf('\n--- Root Locus ---\n');
    [N, D, K0] = getTF('locus');
    m = numel(N) - 1;  n = numel(D) - 1;
    padN = [zeros(1, n-m) N];
    tol = 1e-6;

    z = cleanC(roots(N));
    p = cleanC(roots(D));

    fprintf('\nOpen loop: K * N(S)/D(S)     Characteristic eq: D(S) + K*N(S) = 0,  K >= 0\n');
    fprintf('(on the plot:  x = open-loop poles ,  o = open-loop zeros)\n');
    showPZ(z, p);
    fprintf('\nNumber of branches = %d\n', n);

    % ---------- Asymptotes ----------
    na = n - m;
    angA = [];  sigA = [];
    if na > 0
        angA = (2*(0:na-1) + 1) * 180 / na;
        sigA = (sum(real(p)) - sum(real(z))) / na;
        fprintf('\nAsymptotes (n-m = %d):\n', na);
        fprintf('  Centroid = %.4f\n', sigA);
        fprintf('  Angles   = %s deg\n', mat2str(angA, 5));
    else
        fprintf('\nn = m: no asymptotes (branches end on finite zeros).\n');
    end

    % ---------- Real-axis segments ----------
    realPts = [real(p(abs(imag(p)) < tol)); real(z(abs(imag(z)) < tol))];
    realPts = sort(realPts, 'descend');
    fprintf('\nReal-axis locus segments:\n');
    if isempty(realPts)
        fprintf('  none\n');
    else
        for k = 1:2:numel(realPts)
            if k < numel(realPts)
                fprintf('  [%.4f , %.4f]\n', realPts(k+1), realPts(k));
            else
                fprintf('  (-inf , %.4f]\n', realPts(k));
            end
        end
    end
    onLocus = @(x) mod(sum(realPts > x + 1e-9), 2) == 1;

    % ---------- Break-away / break-in points ----------
    bp = [];  bpK = [];
    B = polySub(conv(polyder(D), N), conv(D, polyder(N)));
    if any(abs(B) > 1e-12)
        B = trimLead(B);
        if numel(B) > 1
            cand = roots(B);
            cand = real(cand(abs(imag(cand)) < tol));
            for c = cand(:)'
                if abs(polyval(N, c)) < 1e-12, continue; end
                Kb = -polyval(D, c) / polyval(N, c);
                if Kb > 1e-9 && onLocus(c)
                    bp(end+1) = c;  bpK(end+1) = Kb; %#ok<AGROW>
                end
            end
        end
    end
    fprintf('\nBreak-away / break-in points:\n');
    if isempty(bp)
        fprintf('  none\n');
    else
        for k = 1:numel(bp)
            fprintf('  S = %.4f   at K = %.4f\n', bp(k), bpK(k));
        end
    end

    % ---------- Angles of departure / arrival ----------
    fprintf('\nAngles of departure (complex poles):\n');
    cp = p(imag(p) > tol);
    if isempty(cp), fprintf('  none (no complex poles)\n'); end
    for k = 1:numel(cp)
        pk = cp(k);
        others = p;  others(find(abs(p - pk) < 1e-9, 1)) = [];
        ang = 180 + sum(angle(pk - z))*180/pi - sum(angle(pk - others))*180/pi;
        ang = mod(ang + 180, 360) - 180;
        fprintf('  pole %s : %.2f deg   (conjugate: %.2f deg)\n', fmtC(pk), ang, -ang);
    end
    fprintf('Angles of arrival (complex zeros):\n');
    cz = z(imag(z) > tol);
    if isempty(cz), fprintf('  none (no complex zeros)\n'); end
    for k = 1:numel(cz)
        zk = cz(k);
        otherz = z;  otherz(find(abs(z - zk) < 1e-9, 1)) = [];
        ang = 180 - sum(angle(zk - otherz))*180/pi + sum(angle(zk - p))*180/pi;
        ang = mod(ang + 180, 360) - 180;
        fprintf('  zero %s : %.2f deg   (conjugate: %.2f deg)\n', fmtC(zk), ang, -ang);
    end

    % ---------- Imaginary-axis crossings ----------
    Dj = D .* (1j).^(n:-1:0);
    Nj = N .* (1j).^(m:-1:0);
    Pi = imag(conv(Dj, conj(Nj)));
    Kall = [];  wall = [];
    if any(abs(Pi) > 1e-12)
        wr = roots(trimLead(Pi));
        if ~isempty(wr)
            wr = real(wr(abs(imag(wr)) < 1e-6 & real(wr) > -1e-9));
            wr = unique(round(max(wr, 0), 8));
            for w = wr(:)'
                Nw = polyval(N, 1j*w);  Dw = polyval(D, 1j*w);
                if abs(Nw) < 1e-9, continue; end
                Kc = -Dw / Nw;
                if abs(imag(Kc)) > 1e-6*max(1, abs(Kc)), continue; end
                if real(Kc) > 1e-9
                    Kall(end+1) = real(Kc);  wall(end+1) = w; %#ok<AGROW>
                end
            end
        end
    end
    [Kcrit, ia] = unique(round(Kall, 8));
    wcr = wall(ia);

    fprintf('\nImaginary-axis crossings:\n');
    if isempty(Kcrit)
        fprintf('  none (poles never cross the jw axis for K > 0)\n');
    else
        for k = 1:numel(Kcrit)
            fprintf('  K = %.4f  at  w = %.4f rad/s\n', Kcrit(k), wcr(k));
        end
    end

    % ---------- Stability condition ----------
    fprintf('\nStability by K range (K > 0):\n');
    edges = [0, Kcrit(:)', inf];
    stableRanges = {};
    for k = 1:numel(edges)-1
        lo = edges(k);  hi = edges(k+1);
        if isinf(hi), pt = max(2*lo, 1); else, pt = (lo + hi)/2; end
        rr = roots(D + pt*padN);
        nR = sum(real(rr) > 1e-7);
        if nR == 0 && all(real(rr) < -1e-7)
            fprintf('  %s : STABLE\n', rangeStr(lo, hi));
            stableRanges{end+1} = rangeStr(lo, hi); %#ok<AGROW>
        else
            fprintf('  %s : UNSTABLE (%d RHP poles)\n', rangeStr(lo, hi), nR);
        end
    end
    if isempty(stableRanges)
        fprintf('Stability condition: system is unstable for every K > 0.\n');
        condTxt = 'unstable for all K>0';
    else
        condTxt = strjoin(stableRanges, ' , ');
        fprintf('Stability condition: stable for  %s\n', condTxt);
    end

    % ---------- Closed loop at the given gain ----------
    clp = cleanC(roots(D + K0*padN));
    fprintf('\nClosed-loop poles at K = %g:\n', K0);
    for k = 1:numel(clp), fprintf('   %s\n', fmtC(clp(k))); end
    if all(real(clp) < -1e-7)
        statusTxt = 'STABLE';
    elseif any(real(clp) > 1e-7)
        statusTxt = 'UNSTABLE';
    else
        statusTxt = 'MARGINALLY STABLE';
    end
    fprintf('At K = %g the closed loop is %s\n', K0, statusTxt);

    % ---------- Plot ----------
    Kmax = 1e4 * max([1, Kcrit, abs(K0), bpK]);
    Kvec = unique([0, logspace(-6, log10(Kmax), 4000), Kcrit, bpK]);
    br = trackLocus(D, padN, Kvec);

    mx = 2 * max([1; abs(p(:)); abs(z(:)); wcr(:); abs(bp(:)); abs(sigA(:))]);

    figure('Name', 'Root Locus', 'Color', 'w', 'Units', 'normalized', ...
           'Position', [0.1 0.08 0.8 0.84]);
    hold on; grid on; box on;
    patch([0 mx mx 0], [-mx -mx mx mx], [1.00 0.90 0.90], 'EdgeColor', 'none');
    patch([-mx 0 0 -mx], [-mx -mx mx mx], [0.91 0.97 0.91], 'EdgeColor', 'none');
    set(gca, 'Layer', 'top');
    plot([-mx mx], [0 0], 'k', 'LineWidth', 0.8);
    plot([0 0], [-mx mx], 'k', 'LineWidth', 0.8);

    colr = lines(max(n, 1));
    for b = 1:n
        plot(real(br(b,:)), imag(br(b,:)), 'Color', colr(b,:), 'LineWidth', 1.8);
    end
    for a = angA
        t = [0 2*mx];
        plot(sigA + t*cosd(a), t*sind(a), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
    end

    hs = [];  nm = {};
    h = plot(real(p), imag(p), 'x', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2.5);
    hs(end+1) = h(1);  nm{end+1} = 'x  Open-loop poles'; %#ok<AGROW>
    if ~isempty(z)
        h = plot(real(z), imag(z), 'o', 'Color', 'k', 'MarkerSize', 9, ...
                 'LineWidth', 2, 'MarkerFaceColor', 'w');
        hs(end+1) = h(1);  nm{end+1} = 'o  Open-loop zeros'; %#ok<AGROW>
    end
    if ~isempty(wcr)
        h = plot(zeros(size(wcr)), wcr, 's', 'MarkerFaceColor', [0.75 0 0.75], ...
                 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
        plot(zeros(size(wcr)), -wcr, 's', 'MarkerFaceColor', [0.75 0 0.75], ...
             'MarkerEdgeColor', 'w', 'MarkerSize', 8);
        hs(end+1) = h(1);  nm{end+1} = 'jw-axis crossing'; %#ok<AGROW>
        for k = 1:numel(wcr)
            text(0.02*mx, wcr(k), sprintf('K=%.3g', Kcrit(k)), 'Color', [0.55 0 0.55], ...
                 'FontSize', 9, 'BackgroundColor', 'w', 'Margin', 1);
            text(0.02*mx, -wcr(k), sprintf('K=%.3g', Kcrit(k)), 'Color', [0.55 0 0.55], ...
                 'FontSize', 9, 'BackgroundColor', 'w', 'Margin', 1);
        end
    end
    if ~isempty(bp)
        h = plot(bp, zeros(size(bp)), '^', 'MarkerFaceColor', [0.1 0.7 0.2], ...
                 'MarkerEdgeColor', 'w', 'MarkerSize', 9);
        hs(end+1) = h(1);  nm{end+1} = 'Break-away / break-in'; %#ok<AGROW>
        for k = 1:numel(bp)
            text(bp(k), 0.05*mx, sprintf('K=%.3g', bpK(k)), 'Color', [0.05 0.45 0.1], ...
                 'FontSize', 9, 'HorizontalAlignment', 'center', ...
                 'BackgroundColor', 'w', 'Margin', 1);
        end
    end
    h = plot(real(clp), imag(clp), 's', 'MarkerFaceColor', [0.85 0.1 0.1], ...
             'MarkerEdgeColor', 'w', 'MarkerSize', 9);
    hs(end+1) = h(1);  nm{end+1} = sprintf('Closed-loop poles (K=%g)', K0); %#ok<AGROW>

    axis equal; xlim([-mx mx]); ylim([-mx mx]);
    xlabel('Real axis (\sigma)'); ylabel('Imaginary axis (j\omega)');
    title(sprintf('Root Locus  |  K = %g -> %s  |  Stable for: %s', K0, statusTxt, condTxt));
    legend(hs, nm, 'Location', 'best');
    text(0.55*mx, -0.93*mx, 'RHP', 'Color', [0.8 0 0], 'FontWeight', 'bold');
    text(-0.85*mx, -0.93*mx, 'LHP', 'Color', [0 0.4 0], 'FontWeight', 'bold');
    hold off;
end

function br = trackLocus(D, padN, Kvec)
    n = numel(D) - 1;
    br = zeros(n, numel(Kvec));
    prev = roots(D + Kvec(1)*padN);
    br(:,1) = prev;
    for i = 2:numel(Kvec)
        r = roots(D + Kvec(i)*padN);
        r = r(:);
        if numel(r) < n, r = [r; 1e6*ones(n-numel(r), 1)]; end
        dist = abs(prev - r.');
        ord = zeros(n, 1);
        for t = 1:n
            [~, idx] = min(dist(:));
            [a, b] = ind2sub([n n], idx);
            ord(a) = b;
            dist(a,:) = inf;  dist(:,b) = inf;
        end
        prev = r(ord);
        br(:,i) = prev;
    end
end

%% ===================================================================
%                              BODE
%% ===================================================================
function runBode()
    fprintf('\n--- Bode Plot ---\n');
    [N, D, K0] = getTF('bode');
    n = numel(D) - 1;  m = numel(N) - 1;
    padN = [zeros(1, n-m) N];

    z = cleanC(roots(N));
    p = cleanC(roots(D));
    fprintf('\nTF = K * N(S)/D(S),  K = %g\n', K0);
    showPZ(z, p);

    % ---------- straight-line (asymptotic) Bode data ----------
    [cw, segSlope, c0, slope0] = asymptoteData(N, D, K0);

    % ---------- frequency range ----------
    if isempty(cw)
        wlo = 0.1;  whi = 100;
    else
        wlo = 10^(floor(log10(min(cw))) - 1);
        whi = 10^(ceil(log10(max(cw))) + 1);
    end
    for it = 1:12
        w = makeGrid(wlo, whi);
        mag = 20*log10(abs(evalTF(N, D, K0, w)));
        changed = false;
        if n > m && mag(end) > -40 && whi < 1e9
            whi = whi * 10;  changed = true;
        end
        if slope0 < 0 && mag(1) < 20 && wlo > 1e-9
            wlo = wlo / 10;  changed = true;
        end
        if ~changed, break; end
    end
    w = makeGrid(wlo, whi);
    H = evalTF(N, D, K0, w);
    mag = 20*log10(abs(H));
    ph = unwrap(angle(H)) * 180/pi;

    % ---------- slopes ----------
    fprintf('\nMagnitude slopes (asymptotic):\n');
    lo = [0, cw];  hi = [cw, inf];
    for k = 1:numel(segSlope)
        if lo(k) == 0 && isinf(hi(k)), rs = 'all w';
        elseif isinf(hi(k)),           rs = sprintf('w > %.4g', lo(k));
        elseif lo(k) == 0,             rs = sprintf('w < %.4g', hi(k));
        else,                          rs = sprintf('%.4g < w < %.4g', lo(k), hi(k));
        end
        fprintf('  %-24s : %s\n', rs, slopeStr(segSlope(k)));
    end
    fprintf('Roll-off at high frequency = %s  (n - m = %d)\n', slopeStr(segSlope(end)), n-m);

    % ---------- crossovers and margins ----------
    [wgc, phgc] = crossings(w, mag, ph);           % GX : |G| = 0 dB
    PM = mod(phgc + 360, 360) - 180;               % phase margin (deg)
    refPh = phgc - PM;                             % the -180 line next to it

    wpc = [];  magpc = [];  PXL = [];              % PX : phase = -180 deg
    for L = -180 + 360*(-3:3)
        [wc, mg] = crossings(w, ph - L, mag);
        wpc = [wpc, wc];  magpc = [magpc, mg];  PXL = [PXL, L*ones(1, numel(wc))]; %#ok<AGROW>
    end
    [wpc, ord] = sort(wpc);  magpc = magpc(ord);  PXL = PXL(ord);
    GM = -magpc;                                    % gain margin (dB)

    fprintf('\nStability margins (unity feedback):\n');
    if isempty(wgc)
        fprintf('  GX (gain crossover)  : none (|G| never equals 0 dB) -> PM undefined\n');
    else
        for k = 1:numel(wgc)
            fprintf('  GX = %.4g rad/s  ->  PM = %.2f deg\n', wgc(k), PM(k));
        end
    end
    if isempty(wpc)
        fprintf('  PX (phase crossover) : none (phase never reaches -180 deg) -> GM = infinite\n');
    else
        for k = 1:numel(wpc)
            fprintf('  PX = %.4g rad/s  ->  GM = %.2f dB\n', wpc(k), GM(k));
        end
    end

    % ---------- UGF and GBW ----------
    if isempty(wgc), UGF = NaN; else, UGF = wgc(1); end
    slopeAtUGF = NaN;
    if ~isnan(UGF)
        slopeAtUGF = segSlope(find(UGF > [0 cw], 1, 'last'));
    end
    GBW = NaN;  gbwNote = '';
    if slope0 == 0                                   % finite, non-zero DC gain
        target = 20*log10(c0) - 3.0103;
        w3 = crossings(w, mag - target, mag);
        if ~isempty(w3)
            GBW = c0 * w3(1);
            gbwNote = sprintf('DC gain %.4g x BW(-3dB) %.4g rad/s', c0, w3(1));
        end
    end
    if isnan(GBW) && ~isnan(UGF)
        GBW = UGF;
        gbwNote = 'equals UGF since |G(UGF)| = 1';
    end

    fprintf('\nUGF / GBW:\n');
    if isnan(UGF)
        fprintf('  UGF : none (|G| never crosses 0 dB)\n');
    else
        fprintf('  UGF = %.4g rad/s (%.4g Hz)   slope at UGF = %s\n', ...
                UGF, UGF/(2*pi), slopeStr(slopeAtUGF));
    end
    if isnan(GBW)
        fprintf('  GBW : not defined\n');
    else
        fprintf('  GBW = %.4g rad/s (%.4g Hz)   [%s]\n', GBW, GBW/(2*pi), gbwNote);
    end

    clp = roots(D + K0*padN);
    if all(real(clp) < -1e-7)
        statusTxt = 'STABLE';
    elseif any(real(clp) > 1e-7)
        statusTxt = 'UNSTABLE';
    else
        statusTxt = 'MARGINALLY STABLE';
    end
    fprintf('Closed loop (unity feedback) is %s\n', statusTxt);

    % ================== PLOT ==================
    colMag = [0.00 0.35 0.75];  colPh  = [0.55 0.15 0.65];
    colAsy = [0.95 0.55 0.10];  colGX  = [0.05 0.60 0.25];
    colPX  = [0.85 0.15 0.15];  colCor = [0.65 0.65 0.65];

    % asymptote vertices
    vw = [wlo, cw, whi];
    vA = zeros(size(vw));
    vA(1) = 20*log10(c0) + slope0*log10(wlo);
    for k = 2:numel(vw)
        vA(k) = vA(k-1) + segSlope(k-1)*(log10(vw(k)) - log10(vw(k-1)));
    end

    % y limits
    fm = mag(isfinite(mag));
    ytop = max([fm, vA]);  ybot = min([fm, vA]);
    ybot = max(ybot, ytop - 300);
    pad = 0.12*(ytop - ybot) + 5;
    yl1 = [ybot - pad, ytop + pad];

    pmin = min([ph, refPh]);  pmax = max([ph, refPh]);
    span = pmax - pmin;
    if span > 900, tick = 180; elseif span > 450, tick = 90; else, tick = 45; end
    yl2 = [tick*floor((pmin - 5)/tick), tick*ceil((pmax + 5)/tick)];

    figure('Name', 'Bode Plot', 'Color', 'w', 'Units', 'normalized', ...
           'Position', [0.08 0.05 0.84 0.85]);

    % ---------------- magnitude ----------------
    ax1 = subplot(2,1,1);
    hold on; grid on; box on;
    set(ax1, 'XScale', 'log');
    xlim([wlo whi]);  ylim(yl1);

    for k = 1:numel(cw)                              % corner frequencies
        plot([cw(k) cw(k)], yl1, ':', 'Color', colCor, 'LineWidth', 1);
        text(cw(k), yl1(2), sprintf(' %.3g', cw(k)), 'Color', [0.4 0.4 0.4], ...
             'FontSize', 8, 'VerticalAlignment', 'top');
    end
    for k = 1:numel(wgc)
        plot([wgc(k) wgc(k)], yl1, ':', 'Color', colGX, 'LineWidth', 1.2);
    end
    for k = 1:numel(wpc)
        plot([wpc(k) wpc(k)], yl1, ':', 'Color', colPX, 'LineWidth', 1.2);
    end
    plot([wlo whi], [0 0], '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1);
    plot(vw, vA, '--', 'Color', colAsy, 'LineWidth', 1.6);          % asymptotes
    plot(w, mag, 'Color', colMag, 'LineWidth', 2.2);                 % actual curve

    sv = [wlo, cw, whi];                                             % slope labels
    for k = 1:numel(segSlope)
        wm = sqrt(sv(k)*sv(k+1));
        ym = vA(k) + segSlope(k)*(log10(wm) - log10(sv(k)));
        if k == numel(segSlope)
            lbl = ['Roll-off ' slopeStr(segSlope(k))];
        else
            lbl = slopeStr(segSlope(k));
        end
        text(wm, ym, [' ' lbl ' '], 'Color', colAsy*0.7, 'FontWeight', 'bold', ...
             'FontSize', 9, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'BackgroundColor', 'w', 'Margin', 1);
    end

    for k = 1:numel(wgc)                                             % GX
        plot(wgc(k), 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', colGX, ...
             'MarkerEdgeColor', 'w', 'LineWidth', 1);
        if numel(wgc) == 1, nmGX = 'GX (UGF)'; else, nmGX = sprintf('GX%d', k); end
        text(wgc(k), 0, sprintf('  %s = %.4g rad/s', nmGX, wgc(k)), 'Color', colGX*0.8, ...
             'FontWeight', 'bold', 'FontSize', 9, 'VerticalAlignment', 'bottom', ...
             'BackgroundColor', 'w', 'Margin', 1);
    end
    for k = 1:numel(wpc)                                             % GM bar at PX
        plot([wpc(k) wpc(k)], [0 magpc(k)], '-', 'Color', colPX, 'LineWidth', 3);
        plot(wpc(k), magpc(k), 's', 'MarkerSize', 8, 'MarkerFaceColor', colPX, ...
             'MarkerEdgeColor', 'w');
        text(wpc(k), magpc(k)/2, {sprintf('  GM = %.2f dB', GM(k)), ...
             sprintf('  (at PX = %.4g)', wpc(k))}, 'Color', colPX*0.85, ...
             'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', 'w', 'Margin', 1);
    end

    if isnan(UGF), s1 = 'UGF = none';
    else, s1 = sprintf('UGF = %.4g rad/s (%.4g Hz)', UGF, UGF/(2*pi)); end
    if isnan(GBW), s2 = 'GBW = n/a';
    else, s2 = sprintf('GBW = %.4g rad/s (%.4g Hz)', GBW, GBW/(2*pi)); end
    info = {s1, s2, ['Roll-off = ' slopeStr(segSlope(end))]};
    text(0.985, 0.95, info, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'top', 'BackgroundColor', [1 1 0.90], ...
         'EdgeColor', [0.7 0.7 0.45], 'FontSize', 9, 'Margin', 5);

    ylabel('Magnitude (dB)');
    title(sprintf('Bode Plot  |  K = %g  |  closed loop (unity feedback): %s', K0, statusTxt));
    legend_h = [plot(NaN, NaN, '-', 'Color', colMag, 'LineWidth', 2.2), ...
                plot(NaN, NaN, '--', 'Color', colAsy, 'LineWidth', 1.6)];
    legend(legend_h, {'Actual |G(j\omega)|', 'Asymptotes (slopes)'}, 'Location', 'southwest');

    % ---------------- phase ----------------
    ax2 = subplot(2,1,2);
    hold on; grid on; box on;
    set(ax2, 'XScale', 'log');
    xlim([wlo whi]);  ylim(yl2);
    set(ax2, 'YTick', yl2(1):tick:yl2(2));

    for k = 1:numel(wgc)
        plot([wgc(k) wgc(k)], yl2, ':', 'Color', colGX, 'LineWidth', 1.2);
    end
    for k = 1:numel(wpc)
        plot([wpc(k) wpc(k)], yl2, ':', 'Color', colPX, 'LineWidth', 1.2);
    end
    for L = -180 + 360*(-3:3)
        if L > yl2(1) && L < yl2(2)
            plot([wlo whi], [L L], '--', 'Color', colPX, 'LineWidth', 1);
        end
    end
    plot(w, ph, 'Color', colPh, 'LineWidth', 2.2);

    for k = 1:numel(wpc)                                             % PX
        plot(wpc(k), PXL(k), 's', 'MarkerSize', 9, 'MarkerFaceColor', colPX, ...
             'MarkerEdgeColor', 'w');
        if numel(wpc) == 1, nmPX = 'PX'; else, nmPX = sprintf('PX%d', k); end
        text(wpc(k), PXL(k), sprintf('  %s = %.4g rad/s', nmPX, wpc(k)), 'Color', colPX*0.85, ...
             'FontWeight', 'bold', 'FontSize', 9, 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'w', 'Margin', 1);
    end
    for k = 1:numel(wgc)                                             % PM bar at GX
        plot([wgc(k) wgc(k)], [refPh(k) phgc(k)], '-', 'Color', colGX, 'LineWidth', 3);
        plot(wgc(k), phgc(k), 'o', 'MarkerSize', 9, 'MarkerFaceColor', colGX, ...
             'MarkerEdgeColor', 'w');
        text(wgc(k), (refPh(k) + phgc(k))/2, {sprintf('  PM = %.2f%s', PM(k), char(176)), ...
             sprintf('  (at GX = %.4g)', wgc(k))}, 'Color', colGX*0.8, ...
             'FontWeight', 'bold', 'FontSize', 9, 'BackgroundColor', 'w', 'Margin', 1);
    end

    ylabel('Phase (deg)');  xlabel('Frequency (rad/s)');
    linkaxes([ax1 ax2], 'x');
    xlim([wlo whi]);
end

function w = makeGrid(wlo, whi)
    w = logspace(log10(wlo), log10(whi), max(2000, ceil(300*log10(whi/wlo))));
end

function H = evalTF(N, D, K0, w)
    H = K0 * polyval(N, 1j*w) ./ polyval(D, 1j*w);
end

function [cw, segSlope, c0, slope0] = asymptoteData(N, D, K0)
    % Corner frequencies, slope of every segment and low-frequency gain.
    zO = trailingZeros(N);  pO = trailingZeros(D);
    Nn = N(1:end-zO);       Dn = D(1:end-pO);
    zr = abs(roots(Nn));    pr = abs(roots(Dn));
    allc = [zr(:); pr(:)];
    chg  = [20*ones(numel(zr), 1); -20*ones(numel(pr), 1)];
    [allc, ix] = sort(allc);  chg = chg(ix);

    cw = [];  dchg = [];
    for k = 1:numel(allc)
        if ~isempty(cw) && abs(allc(k) - cw(end)) <= 1e-6*cw(end)
            dchg(end) = dchg(end) + chg(k);
        else
            cw(end+1) = allc(k);  dchg(end+1) = chg(k); %#ok<AGROW>
        end
    end
    slope0   = 20*(zO - pO);
    segSlope = slope0 + [0, cumsum(dchg)];
    c0       = abs(K0 * Nn(end) / Dn(end));      % |G| ~ c0 * w^(slope0/20) at low w
end

function k = trailingZeros(v)
    k = 0;
    while k < numel(v) - 1 && abs(v(end-k)) < 1e-12
        k = k + 1;
    end
end

function s = slopeStr(v)
    v = round(v);
    if v == 0, s = '0 dB/dec'; else, s = sprintf('%+d dB/dec', v); end
end

function [wc, ex] = crossings(w, y, extra)
    i = find(y(1:end-1) .* y(2:end) < 0);
    t = y(i) ./ (y(i) - y(i+1));
    wc = 10.^(log10(w(i)) + t .* (log10(w(i+1)) - log10(w(i))));
    ex = extra(i) + t .* (extra(i+1) - extra(i));
end

%% ===================================================================
%                     SAFE INPUT HELPERS
%% ===================================================================
function v = readNum(msg)
    % Reads one real finite number. Accepts things like 2, -3.5, 1e3, 1/3, 2*pi.
    % Never crashes on bad input: it just asks again.
    while true
        str = strtrim(input(msg, 's'));
        if isempty(str)
            disp('Empty input. Please enter a number.');
            continue;
        end
        v = str2double(str);
        if isnan(v) && ~isempty(regexp(str, '^[0-9eEpi\.\+\-\*/\^\(\)\s]+$', 'once'))
            try
                v = eval(str);
            catch
                v = NaN;
            end
        end
        if isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v)
            return;
        end
        disp('Invalid input. Enter a single real number (e.g. 2, -3.5, 1/3, 2*pi).');
    end
end

function n = readInt(msg, lo, hi)
    while true
        n = readNum(msg);
        if n == floor(n) && n >= lo && n <= hi, return; end
        fprintf('Please enter an integer between %d and %d.\n', lo, hi);
    end
end

function c = readChoice(msg, allowed)
    while true
        c = readNum(msg);
        if ismember(c, allowed), return; end
        fprintf('Please choose one of: %s\n', mat2str(allowed));
    end
end

function [num, den, K0] = getTF(kind)
    % kind = 'locus' (needs K >= 0) or 'bode' (needs K ~= 0)
    m = readInt('Enter order of the NUMERATOR m (0 for a constant) = ', 0, 30);
    n = readInt('Enter order of the DENOMINATOR n = ', 1, 30);
    while m > n
        disp('Improper TF: numerator order must not exceed denominator order.');
        m = readInt('Enter order of the NUMERATOR m (0 for a constant) = ', 0, 30);
        n = readInt('Enter order of the DENOMINATOR n = ', 1, 30);
    end

    fprintf('\nTF = K * N(S)/D(S)   (enter coefficients from the highest power)\n');
    while true
        num = zeros(1, m+1);  den = zeros(1, n+1);
        for k = 1:m+1
            num(k) = readNum(sprintf('Numerator coefficient of S^%d = ', m-k+1));
        end
        for k = 1:n+1
            den(k) = readNum(sprintf('Denominator coefficient of S^%d = ', n-k+1));
        end
        if ~any(num ~= 0)
            disp('The numerator is all zeros. Please enter the coefficients again.');
            continue;
        end
        if ~any(den ~= 0)
            disp('The denominator is all zeros. Please enter the coefficients again.');
            continue;
        end
        num = trimLead(num);
        den = trimLead(den);
        if numel(den) < 2
            disp('The denominator must be of order 1 or higher. Please enter again.');
            continue;
        end
        if numel(num) > numel(den)
            disp('Improper TF: numerator order exceeds denominator order. Please enter again.');
            continue;
        end
        break;
    end

    while true
        K0 = readNum('Enter the gain K = ');
        if strcmp(kind, 'locus') && K0 < 0
            disp('For the root locus, K must be >= 0.');
        elseif strcmp(kind, 'bode') && K0 == 0
            disp('For a Bode plot, K must be non-zero.');
        else
            break;
        end
    end
end

%% ===================================================================
%                          SHARED HELPERS
%% ===================================================================
function v = trimLead(v)
    idx = find(abs(v) > 1e-12, 1);
    if isempty(idx), error('Polynomial cannot be all zeros.'); end
    v = v(idx:end);
end

function c = polySub(a, b)
    L = max(numel(a), numel(b));
    c = [zeros(1, L-numel(a)) a] - [zeros(1, L-numel(b)) b];
end

function x = cleanC(x)
    x = x(:);
    x(abs(imag(x)) < 1e-9) = real(x(abs(imag(x)) < 1e-9));
end

function s = fmtC(x)
    s = sprintf('%.4f %+.4fj', real(x), imag(x));
end

function showPZ(z, p)
    fprintf('\nPoles (%d):\n', numel(p));
    for k = 1:numel(p), fprintf('   p%d = %s\n', k, fmtC(p(k))); end
    fprintf('Zeros (%d):\n', numel(z));
    if isempty(z), fprintf('   none (finite)\n'); end
    for k = 1:numel(z), fprintf('   z%d = %s\n', k, fmtC(z(k))); end
end

function c = countChanges(sg)
    c = sum(sg(1:end-1) .* sg(2:end) < 0);
end

function s = rangeStr(lo, hi)
    if isinf(lo) && isinf(hi), s = 'all K';
    elseif isinf(lo),          s = sprintf('K < %.6g', hi);
    elseif isinf(hi),          s = sprintf('K > %.6g', lo);
    else,                      s = sprintf('%.6g < K < %.6g', lo, hi);
    end
end

%% ===================================================================
%              ROUTH ENGINE (with all special cases)
%% ===================================================================
function A = analyzePoly(S)
    % Runs the numeric Routh test and counts RHP / jw-axis / LHP poles.
    S = trimLead(S);
    n = numel(S) - 1;
    if n < 1, error('Polynomial must be of order 1 or higher.'); end
    [R, special, info] = routhNumeric(S);
    sg = sign(R(:,1));
    numRHP = countChanges(sg);
    numImag = 0;
    for i = info.zeroRows
        degAux = n - i + 2;                                  % degree of the auxiliary polynomial
        numImag = numImag + (degAux - 2*countChanges(sg(i-1:end)));
    end
    numImag = max(numImag, 0);
    A = struct('S', S, 'n', n, 'R', R, 'special', special, 'info', info, ...
               'numRHP', numRHP, 'numImag', numImag, 'numLHP', n - numRHP - numImag);
end

function [R, special, info] = routhNumeric(S)
    % Special case 1: first element of a row is 0 (row not all zero) -> replace by eps.
    % Special case 2: whole row is 0 -> replace by the derivative of the auxiliary polynomial.
    n = numel(S) - 1;
    cols = ceil((n+1)/2) + 1;
    scale = max(abs(S));
    epsVal = 1e-6 * scale;
    special = false;
    info = struct('zeroRows', [], 'epsRows', [], 'aux', {{}}, 'epsVal', epsVal);

    R = zeros(n+1, cols);
    r1 = S(1:2:end);  R(1, 1:numel(r1)) = r1;
    r2 = S(2:2:end);  R(2, 1:numel(r2)) = r2;

    for i = 2:n+1
        if i >= 3
            for j = 1:cols-1
                R(i,j) = (R(i-1,1)*R(i-2,j+1) - R(i-2,1)*R(i-1,j+1)) / R(i-1,1);
            end
        end
        tolz = 1e-9 * max(max(abs(R(max(i-2,1):i-1, :))));   % tolerance relative to the rows above

        if all(abs(R(i,:)) <= tolz)                  % ---- zero row ----
            special = true;
            d = n - i + 2;                           % degree of auxiliary polynomial A(S)
            aux = zeros(1, d+1);
            for j = 1:cols
                idx = 1 + 2*(j-1);
                if idx <= d+1, aux(idx) = R(i-1,j); end
            end
            info.zeroRows(end+1) = i;
            info.aux{end+1} = aux;
            for j = 1:cols                           % dA/dS coefficients
                R(i,j) = (d - 2*(j-1)) * R(i-1,j);
            end
        end

        if abs(R(i,1)) <= tolz && i < n+1            % ---- zero first element ----
            special = true;
            R(i,1) = epsVal;
            info.epsRows(end+1) = i;
        end
    end
end

function printRouthTable(R, n, info)
    cols = size(R, 2);
    fprintf('\nRouth Table:\n');
    for i = 1:n+1
        fprintf('S^%-2d |', n-i+1);
        fprintf(' %12.4g', R(i, 1:cols-1));
        if ismember(i, info.epsRows)
            fprintf('   <- first element was 0, replaced by eps = %.2g', info.epsVal);
        end
        if ismember(i, info.zeroRows)
            fprintf('   <- zero row, replaced by dA/dS');
        end
        fprintf('\n');
    end
end

function printAux(info, n)
    for k = 1:numel(info.aux)
        a = info.aux{k};
        i = info.zeroRows(k);
        fprintf('\nZero row at S^%d -> auxiliary polynomial A(S) (degree %d) coefficients: %s\n', ...
                n-i+1, numel(a)-1, mat2str(a, 5));
        rr = roots(a);
        fprintf('Roots of A(S):\n');
        for q = 1:numel(rr), fprintf('   %s\n', fmtC(rr(q))); end
    end
end

function tf = hasRepeatedImag(r)
    % true if two (or more) poles sit on top of each other on the jw axis
    ir = r(abs(real(r)) <= 1e-6);
    tf = false;
    for a = 1:numel(ir)-1
        for b = a+1:numel(ir)
            if abs(ir(a) - ir(b)) < 1e-3*max(1, abs(ir(a)))
                tf = true;  return;
            end
        end
    end
end

function plotPoles(r, resultTxt, numRHP)
    tol = 1e-6;
    figure('Name', 'S-domain', 'Color', 'w'); hold on; grid on; box on;
    m = max([abs(real(r)); abs(imag(r)); 1]) * 1.3;
    patch([0 m m 0], [-m -m m m], [1.00 0.90 0.90], 'EdgeColor', 'none');
    patch([-m 0 0 -m], [-m -m m m], [0.91 0.97 0.91], 'EdgeColor', 'none');
    set(gca, 'Layer', 'top');
    plot([-m m], [0 0], 'k', 'LineWidth', 1);
    plot([0 0], [-m m], 'k', 'LineWidth', 1);

    groups = {real(r) < -tol,       [0.00 0.45 0.20], 'Stable poles (LHP)'; ...
              abs(real(r)) <= tol,  [0.95 0.55 0.05], 'Poles on j\omega axis'; ...
              real(r) > tol,        [0.85 0.10 0.10], 'Unstable poles (RHP)'};
    hs = [];  nm = {};
    for g = 1:3
        idx = groups{g,1};
        if any(idx)
            h = plot(real(r(idx)), imag(r(idx)), 'x', 'Color', groups{g,2}, ...
                     'MarkerSize', 12, 'LineWidth', 2.5);
            hs(end+1) = h(1);  nm{end+1} = ['x  ' groups{g,3}]; %#ok<AGROW>
        end
    end
    for k = 1:numel(r)
        text(real(r(k)) + 0.03*m, imag(r(k)) + 0.05*m, ...
            sprintf('%.2f%+.2fj', real(r(k)), imag(r(k))), 'FontSize', 9, ...
            'BackgroundColor', 'w', 'Margin', 1);
    end
    axis equal; xlim([-m m]); ylim([-m m]);
    xlabel('Real axis (\sigma)'); ylabel('Imaginary axis (j\omega)');
    title(sprintf('S-domain poles  |  %s  |  RHP poles = %d', resultTxt, numRHP));
    legend(hs, nm, 'Location', 'best');
    text(0.55*m, -0.9*m, 'RHP', 'Color', [0.8 0 0], 'FontWeight', 'bold');
    text(-0.8*m, -0.9*m, 'LHP', 'Color', [0 0.4 0], 'FontWeight', 'bold');
    hold off;
end
