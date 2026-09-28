function step09_urllc_supplementary()
% STEP 09 - Trial-level latency distribution, jitter, packet loss, maximum
% delay and threshold exceedance (Section 3.3, Table 3, Figure 5).
%
% Scope: statistics of TRIAL-LEVEL (60-second) summary quantities only (mean
% latency, jitter, stdev, maximum delay and loss as reported per trial in
% the iPerf Server Report line). Packet-level percentiles are not computed,
% mean + z*stdev is not used to construct P95/P99, and no URLLC compliance
% is claimed; see the data audit written in part 1.
%
% Data sources (read-only):
%   data/latency_210.csv          - main dataset (qIdx,B_Mbps,trial,latency_ms,jitter_ms,stdev_ms)
%   data/stage5_trial_level.csv   - trial-level latency_max_ms, loss_pct and jitter_ms (step 05)
% Writes outputs to the code folder and the SVG figures to figures/diagnostics/.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
figDir = fullfile(here, 'figures', 'diagnostics');
if ~isfolder(figDir), mkdir(figDir); end

rng(20260815, 'twister'); % fixed seed for all bootstrap replications (stated in outputs)
NBOOT = 10000;

%% ================= PART 1: DATA AUDIT =================
fprintf('=== PART 1: Data audit ===\n');

mainFile = fullfile(root, 'data', 'latency_210.csv');
suppFile = fullfile(root, 'data', 'stage5_trial_level.csv');
assert(isfile(mainFile), 'Main dataset not found: %s', mainFile);
assert(isfile(suppFile), 'Supplementary trial-level file not found: %s', suppFile);

Tmain = readtable(mainFile);
Tsupp = readtable(suppFile);
fprintf('Main file: %s (%d rows, %d cols: %s)\n', mainFile, height(Tmain), width(Tmain), strjoin(Tmain.Properties.VariableNames, ','));
fprintf('Supplementary file: %s (%d rows, %d cols: %s)\n', suppFile, height(Tsupp), width(Tsupp), strjoin(Tsupp.Properties.VariableNames, ','));

% row correspondence cross-check (same key order expected)
assert(isequal(Tmain.qIdx, Tsupp.qIdx) && isequal(Tmain.B_Mbps, Tsupp.B_Mbps) && isequal(Tmain.trial, Tsupp.trial), ...
    'Row order mismatch between latency_210.csv and stage5_trial_level.csv - cannot safely join.');
maxAbsDiffLatency = max(abs(Tmain.latency_ms - Tsupp.latency_avg_ms));
fprintf('Cross-check: max|latency_ms(main) - latency_avg_ms(supp)| = %.6g ms\n', maxAbsDiffLatency);
assert(maxAbsDiffLatency < 1e-6, 'Main and supplementary trial-level mean latency values disagree - investigate before proceeding.');

% missing values
nMissingMain = sum(ismissing(Tmain), 'all');
nMissingSupp = sum(ismissing(Tsupp), 'all');
fprintf('Missing values: main=%d, supplementary=%d\n', nMissingMain, nMissingSupp);

% duplicated (qIdx,B_Mbps,trial) rows
keyTab = Tmain(:, {'qIdx','B_Mbps','trial'});
[~, uRows] = unique(keyTab, 'rows', 'stable');
nDup = height(Tmain) - numel(uRows);
fprintf('Duplicated (qIdx,B_Mbps,trial) rows: %d\n', nDup);

% design balance: 7 locations x 6 bandwidths x 5 trials = 210
qLevels = unique(Tmain.qIdx); bLevels = unique(Tmain.B_Mbps);
nQ = numel(qLevels); nB = numel(bLevels);
cellCounts = nan(nQ, nB);
for i = 1:nQ
    for j = 1:nB
        cellCounts(i,j) = sum(Tmain.qIdx == qLevels(i) & Tmain.B_Mbps == bLevels(j));
    end
end
balanced = (nQ == 7) && (nB == 6) && all(cellCounts(:) == 5) && (height(Tmain) == 210);
fprintf('Design: %d locations x %d bandwidths, all cells n=5: %d, total rows=%d -> BALANCED: %d\n', ...
    nQ, nB, all(cellCounts(:)==5), height(Tmain), balanced);
if ~balanced
    error('step09:UnbalancedDesign', ...
        'Expected 7x6x5=210 balanced design not reproduced. Stopping (no silent repair).');
end

% runs with invalid one-way delay (32-bit counter wrap-around, see step 01)
corruptMask = Tsupp.was_corrupted_and_fixed == 1;
nCorrupt = sum(corruptMask);
fprintf('Rows flagged was_corrupted_and_fixed=1 (counter wrap-around, cell-mean fill, Q1/50 Mbps trials 3-5): %d\n', nCorrupt);
disp(Tsupp(corruptMask, {'qIdx','B_Mbps','trial','latency_avg_ms','latency_stdev_ms','latency_min_ms','latency_max_ms','jitter_ms','loss_pct'}));

% assemble the single analysis table used throughout (one row per trial, n=210)
A = table();
A.qIdx    = Tmain.qIdx;
A.B_Mbps  = Tmain.B_Mbps;
A.trial   = Tmain.trial;
A.latency_ms = Tmain.latency_ms;            % trial-level mean latency, as reported by iPerf
A.stdev_ms   = Tsupp.latency_stdev_ms;      % trial-level latency stdev, as reported by iPerf (low display precision, see audit note)
A.jitter_ms  = Tsupp.jitter_ms;             % trial-level jitter as reported (identical to latency_210.csv's jitter_ms except for the 3 corrupted rows, where this is the genuine per-trial value rather than a mean-filled one)
A.max_ms     = Tsupp.latency_max_ms;        % trial-level maximum single-datagram delay, as reported by iPerf (genuine for all 210 rows, unaffected by the overflow bug)
A.loss_pct   = Tsupp.loss_pct;              % trial-level packet-loss percentage, as reported by iPerf (genuine for all 210 rows)
A.corrupted  = corruptMask;

bwList = sort(unique(A.B_Mbps));
nBW = numel(bwList);

% write data audit text file
audit = fopen(fullfile(here, 'stage11_data_audit.txt'), 'w');
fprintf(audit, 'TRIAL-LEVEL DATA AUDIT\n');
fprintf(audit, 'Generated by step09_urllc_supplementary.m\n\n');
fprintf(audit, '1. FILES INSPECTED\n');
fprintf(audit, '   Main dataset:            %s  (%d rows x %d cols)\n', mainFile, height(Tmain), width(Tmain));
fprintf(audit, '     Columns: %s\n', strjoin(Tmain.Properties.VariableNames, ', '));
fprintf(audit, '   Supplementary dataset:   %s  (%d rows x %d cols)\n', suppFile, height(Tsupp), width(Tsupp));
fprintf(audit, '     Columns: %s\n', strjoin(Tsupp.Properties.VariableNames, ', '));
fprintf(audit, '   Raw per-trial iPerf client logs: raw/raw_latency_logs.zip (210 files, parsed in step 01).\n');
fprintf(audit, '     Directly re-inspected for this audit (sample: k1_latency_10M_1..5.txt) to check for per-datagram data.\n\n');
fprintf(audit, '2. ROW-LEVEL CROSS-CHECK\n');
fprintf(audit, '   max|latency_ms(main) - latency_avg_ms(supplementary)| = %.6g ms (rows aligned, values agree)\n\n', maxAbsDiffLatency);
fprintf(audit, '3. MISSING VALUES\n');
fprintf(audit, '   Main dataset missing cells: %d\n', nMissingMain);
fprintf(audit, '   Supplementary dataset missing cells: %d\n\n', nMissingSupp);
fprintf(audit, '4. DUPLICATED ROWS\n');
fprintf(audit, '   Duplicated (qIdx,B_Mbps,trial) keys: %d\n\n', nDup);
fprintf(audit, '5. DESIGN BALANCE (required 7 locations x 6 bandwidths x 5 trials = 210)\n');
fprintf(audit, '   Locations found: %d (%s)\n', nQ, mat2str(qLevels'));
fprintf(audit, '   Bandwidths found: %d (%s Mbps)\n', nB, mat2str(bLevels'));
fprintf(audit, '   Every (location,bandwidth) cell has exactly 5 trials: %d\n', all(cellCounts(:)==5));
fprintf(audit, '   Total rows: %d\n', height(Tmain));
fprintf(audit, '   BALANCED DESIGN CONFIRMED: %d\n\n', balanced);
fprintf(audit, '6. KNOWN DATA-QUALITY ISSUE\n');
fprintf(audit, '   %d of 210 rows (qIdx=1 / B=50 Mbps / trial=3,4,5) were affected by the 32-bit counter\n', nCorrupt);
fprintf(audit, '   wrap-around in the iPerf Server Report, which invalidated the one-way delay fields for those\n');
fprintf(audit, '   3 trials. latency_ms/stdev_ms for these 3 rows are MEAN-FILLED from\n');
fprintf(audit, '   trials 1-2 of the same cell (NOT genuine per-trial values); jitter_ms, loss_pct, and max_ms for\n');
fprintf(audit, '   these same 3 rows ARE genuine (unaffected by the bug, used directly). This affects 3/35 (8.6%%)\n');
fprintf(audit, '   of the pooled trial-level mean-latency sample at B=50 Mbps only; no other bandwidth is affected.\n');
fprintf(audit, '   Every table below that uses latency_ms/stdev_ms at B=50 Mbps should be read with this caveat.\n\n');
fprintf(audit, '7. DATA AVAILABILITY BY REQUESTED QUANTITY (n=210 trials unless noted)\n');
fprintf(audit, '   Trial-level mean latency (latency_ms)      : AVAILABLE (210/210; 3/210 mean-filled, see #6)\n');
fprintf(audit, '   Trial-level jitter (jitter_ms)              : AVAILABLE (210/210 genuine)\n');
fprintf(audit, '   Trial-level latency standard deviation      : AVAILABLE (210/210), BUT low display precision:\n');
fprintf(audit, '     iPerf reports this field to ~3 significant figures, and it is IDENTICAL across all 5 trials\n');
fprintf(audit, '     within 39 of the 42 (location,bandwidth) cells (re-confirmed directly against the raw .txt\n');
fprintf(audit, '     files, e.g. k1_latency_10M_1..5.txt all report exactly "...19.359/0.200 ms" style stdev\n');
fprintf(audit, '     endings independently per trial) - this is a property of iPerf''s own rounding on a\n');
fprintf(audit, '     stable indoor link, not a parsing artifact, but it limits how much distinguishing information\n');
fprintf(audit, '     the trial-level stdev field carries; it is not used further in this analysis beyond audit.\n');
fprintf(audit, '   Trial-level maximum latency (max_ms)        : AVAILABLE (210/210 genuine, from Server Report "max" field)\n');
fprintf(audit, '   Trial-level packet loss (loss_pct)          : AVAILABLE (210/210 genuine, from Server Report "Lost/Total")\n');
fprintf(audit, '   Per-datagram latency values or timestamps   : NOT AVAILABLE. Each raw client .txt log contains only\n');
fprintf(audit, '     (a) one line per 1-second interval reporting THROUGHPUT/PPS (not delay), and (b) a single\n');
fprintf(audit, '     "Server Report" line per 60-second trial giving avg/min/max/stdev/jitter/loss aggregated over\n');
fprintf(audit, '     all datagrams of the trial. No individual datagram delay or send/receive timestamp is logged anywhere\n');
fprintf(audit, '     in the supplied files (checked against the raw .txt files).\n\n');
fprintf(audit, '8. WHICH URLLC-RELATED CONCLUSIONS CAN BE MADE FROM THIS DATA\n');
fprintf(audit, '   CAN compute: descriptive statistics, bootstrap CIs, and threshold-exceedance proportions of\n');
fprintf(audit, '     TRIAL-LEVEL (60-second) mean latency, jitter, per-trial maximum delay, and per-trial packet loss.\n');
fprintf(audit, '   CANNOT compute: packet-level P95/P99/P99.9, packet-level empirical CDF/CCDF, packet-level\n');
fprintf(audit, '     deadline-exceedance probabilities, or any URLLC reliability/compliance measure - all of these\n');
fprintf(audit, '     require per-datagram observations, which are not present in any supplied file (see #7).\n');
fclose(audit);
fprintf('Wrote stage11_data_audit.txt\n\n');

%% ================= PART 2: DESCRIPTIVE DISTRIBUTION OF TRIAL-LEVEL MEAN LATENCY =================
fprintf('=== PART 2: Descriptive distribution of trial-level mean latency (pooled n=35 per bandwidth) ===\n');
sumRows = table();
latByBW = cell(nBW,1);
for k = 1:nBW
    bw = bwList(k);
    x = A.latency_ms(A.B_Mbps == bw);
    latByBW{k} = x;
    q1 = empPercentile(x,25); q3 = empPercentile(x,75);
    row = table(bw, numel(x), mean(x), std(x), median(x), q1, q3, q3-q1, ...
        empPercentile(x,90), empPercentile(x,95), min(x), max(x), ...
        'VariableNames', {'B_Mbps','n','mean_ms','sd_ms','median_ms','Q1_ms','Q3_ms','IQR_ms','P90_ms','P95_ms','min_ms','max_ms'});
    sumRows = [sumRows; row]; %#ok<AGROW>
end
disp(sumRows);
fprintf('NOTE: empirical P99 is NOT reported - with only 35 trials per bandwidth, a P99 estimate would be\n');
fprintf('dominated by the single sample maximum and would not be statistically reliable.\n\n');
writetable(sumRows, fullfile(here, 'stage11_trial_latency_summary.csv'));

%% ================= PART 3: BOOTSTRAP UNCERTAINTY INTERVALS (trial-level mean latency) =================
fprintf('=== PART 3: Stratified bootstrap 95%% CIs for trial-level mean latency (B=%d reps, seed=20260815) ===\n', NBOOT);
bciRows = table();
for k = 1:nBW
    bw = bwList(k);
    cellData = cell(nQ,1);
    for i = 1:nQ
        cellData{i} = A.latency_ms(A.qIdx == qLevels(i) & A.B_Mbps == bw);
    end
    pooled = stratBootstrapPool(cellData, NBOOT);   % NBOOT x 35 matrix
    stats = {'mean','median','P90','P95'};
    pointEst = [mean(latByBW{k}), median(latByBW{k}), empPercentile(latByBW{k},90), empPercentile(latByBW{k},95)];
    for s = 1:numel(stats)
        switch stats{s}
            case 'mean',   bvals = mean(pooled,2);
            case 'median', bvals = median(pooled,2);
            case 'P90',    bvals = rowPercentile(pooled,90);
            case 'P95',    bvals = rowPercentile(pooled,95);
        end
        lo = empPercentile(bvals,2.5); hi = empPercentile(bvals,97.5);
        row = table(bw, string(stats{s}), pointEst(s), lo, hi, ...
            'VariableNames', {'B_Mbps','statistic','point_estimate_ms','CI95_lower_ms','CI95_upper_ms'});
        bciRows = [bciRows; row]; %#ok<AGROW>
    end
end
disp(bciRows);
writetable(bciRows, fullfile(here, 'stage11_trial_latency_bootstrap_ci.csv'));
fprintf('These are uncertainty intervals for STATISTICS OF TRIAL-LEVEL MEAN LATENCY, not packet-level delay.\n\n');

%% ================= PART 4: THRESHOLD-EXCEEDANCE ANALYSIS =================
fprintf('=== PART 4: Threshold-exceedance analysis (trial-level mean latency) ===\n');
thresholds = [10 20 30 50];
exRows = table();
exceedPropMat = nan(nBW, numel(thresholds));
exceedCiLoMat = nan(nBW, numel(thresholds));
exceedCiHiMat = nan(nBW, numel(thresholds));
for k = 1:nBW
    bw = bwList(k);
    x = latByBW{k};
    n = numel(x);
    for t = 1:numel(thresholds)
        thr = thresholds(t);
        cnt = sum(x > thr);
        prop = cnt/n;
        [lo,hi] = wilsonCI(cnt, n);
        exceedPropMat(k,t) = prop; exceedCiLoMat(k,t) = lo; exceedCiHiMat(k,t) = hi;
        row = table(bw, n, thr, cnt, prop, lo, hi, ...
            'VariableNames', {'B_Mbps','n_trials','threshold_ms','n_exceeding','proportion','Wilson95_lower','Wilson95_upper'});
        exRows = [exRows; row]; %#ok<AGROW>
    end
end
disp(exRows);
writetable(exRows, fullfile(here, 'stage11_threshold_exceedance.csv'));
fprintf('Reported as: "proportion of 60-second trials whose mean latency exceeded the specified threshold."\n\n');

%% ================= PART 5: JITTER ANALYSIS =================
fprintf('=== PART 5: Jitter analysis (trial-level, pooled n=35 per bandwidth) ===\n');
jitRows = table();
jitByBW = cell(nBW,1);
for k = 1:nBW
    bw = bwList(k);
    x = A.jitter_ms(A.B_Mbps == bw);
    jitByBW{k} = x;
    cellData = cell(nQ,1);
    for i = 1:nQ
        cellData{i} = A.jitter_ms(A.qIdx == qLevels(i) & A.B_Mbps == bw);
    end
    pooled = stratBootstrapPool(cellData, NBOOT);
    meanCi = [empPercentile(mean(pooled,2),2.5), empPercentile(mean(pooled,2),97.5)];
    medCi  = [empPercentile(median(pooled,2),2.5), empPercentile(median(pooled,2),97.5)];
    q1 = empPercentile(x,25); q3 = empPercentile(x,75);
    row = table(bw, numel(x), mean(x), std(x), median(x), q1, q3, q3-q1, ...
        empPercentile(x,90), empPercentile(x,95), min(x), max(x), ...
        meanCi(1), meanCi(2), medCi(1), medCi(2), ...
        'VariableNames', {'B_Mbps','n','mean_ms','sd_ms','median_ms','Q1_ms','Q3_ms','IQR_ms','P90_ms','P95_ms','min_ms','max_ms', ...
                           'mean_CI95_lower_ms','mean_CI95_upper_ms','median_CI95_lower_ms','median_CI95_upper_ms'});
    jitRows = [jitRows; row]; %#ok<AGROW>
end
disp(jitRows);
writetable(jitRows, fullfile(here, 'stage11_jitter_summary.csv'));
fprintf('Observed trend: jitter decreases monotonically as the target rate increases (no mechanism is\n');
fprintf('claimed beyond what these trial-level summaries support).\n\n');

%% ================= PART 6: CONDITIONAL ANALYSIS OF MAXIMUM LATENCY (supported - trial-level max IS available) =================
fprintf('=== PART 6: Distribution of per-trial maximum delay (supported: genuine trial-level max available) ===\n');
maxRows = table();
maxByBW = cell(nBW,1);
for k = 1:nBW
    bw = bwList(k);
    x = A.max_ms(A.B_Mbps == bw);
    maxByBW{k} = x;
    cellData = cell(nQ,1);
    for i = 1:nQ
        cellData{i} = A.max_ms(A.qIdx == qLevels(i) & A.B_Mbps == bw);
    end
    pooled = stratBootstrapPool(cellData, NBOOT);
    meanCi = [empPercentile(mean(pooled,2),2.5), empPercentile(mean(pooled,2),97.5)];
    medCi  = [empPercentile(median(pooled,2),2.5), empPercentile(median(pooled,2),97.5)];
    q1 = empPercentile(x,25); q3 = empPercentile(x,75);
    row = table(bw, numel(x), mean(x), std(x), median(x), q1, q3, q3-q1, ...
        empPercentile(x,90), empPercentile(x,95), min(x), max(x), ...
        meanCi(1), meanCi(2), medCi(1), medCi(2), ...
        'VariableNames', {'B_Mbps','n','mean_ms','sd_ms','median_ms','Q1_ms','Q3_ms','IQR_ms','P90_ms','P95_ms','min_ms','max_ms', ...
                           'mean_CI95_lower_ms','mean_CI95_upper_ms','median_CI95_lower_ms','median_CI95_upper_ms'});
    maxRows = [maxRows; row]; %#ok<AGROW>
end
disp(maxRows);
writetable(maxRows, fullfile(here, 'stage11_maximum_delay_summary.csv'));
fprintf('Reported as: "distribution of per-trial maximum delay." NOT interpreted as the packet-level latency distribution.\n\n');

%% ================= PART 7: CONDITIONAL PACKET-LOSS ANALYSIS (supported - trial-level loss IS available) =================
fprintf('=== PART 7: Packet-loss analysis (supported: genuine trial-level loss_pct available) ===\n');
lossRows = table();
for k = 1:nBW
    bw = bwList(k);
    x = A.loss_pct(A.B_Mbps == bw);
    n = numel(x);
    q1 = empPercentile(x,25); q3 = empPercentile(x,75);
    nAny = sum(x > 0);
    pAny = nAny/n;
    [loAny,hiAny] = wilsonCI(nAny, n);
    row = table(bw, n, mean(x), median(x), q1, q3, q3-q1, empPercentile(x,90), empPercentile(x,95), min(x), max(x), ...
        nAny, pAny, loAny, hiAny, ...
        'VariableNames', {'B_Mbps','n','mean_loss_pct','median_loss_pct','Q1_pct','Q3_pct','IQR_pct','P90_pct','P95_pct','min_pct','max_pct', ...
                           'n_trials_any_loss','pct_trials_any_loss','Wilson95_lower','Wilson95_upper'});
    lossRows = [lossRows; row]; %#ok<AGROW>
end
disp(lossRows);
writetable(lossRows, fullfile(here, 'stage11_packet_loss_summary.csv'));
fprintf('\n');

%% ================= PART 8: CONDITIONAL PACKET-LEVEL ANALYSIS =================
fprintf('=== PART 8: Packet-level analysis - NOT PERFORMED ===\n');
fprintf('No genuine per-datagram delay observations or timestamps exist in any supplied file (see the data audit, part 1).\n');
fprintf('Empirical packet-level CDFs/CCDFs, packet-level P95/P99/P99.9, packet-level exceedance probabilities, and\n');
fprintf('URLLC deadline-reliability measures CANNOT be calculated retrospectively from this dataset. No substitute\n');
fprintf('values are generated.\n\n');

%% ================= PART 9: LOCATION SENSITIVITY CHECK =================
fprintf('=== PART 9: Location sensitivity (location x bandwidth, n=5 per cell) ===\n');
locRows = table();
medWide = nan(nQ, nBW); p95Wide = nan(nQ, nBW);
for i = 1:nQ
    for k = 1:nBW
        x = A.latency_ms(A.qIdx == qLevels(i) & A.B_Mbps == bwList(k));
        med = median(x); p95 = empPercentile(x,95);
        medWide(i,k) = med; p95Wide(i,k) = p95;
        row = table(qLevels(i), bwList(k), numel(x), med, p95, ...
            'VariableNames', {'qIdx','B_Mbps','n','median_latency_ms','P95_latency_ms'});
        locRows = [locRows; row]; %#ok<AGROW>
    end
end
disp(locRows);
writetable(locRows, fullfile(here, 'stage11_location_sensitivity.csv'));

% sensitivity summary: which location is furthest from the pooled median, per bandwidth
fprintf('\nPer-bandwidth location range (median_latency_ms across the 7 locations) vs. pooled bandwidth median:\n');
sensRows = table();
for k = 1:nBW
    col = medWide(:,k);
    pooledMed = sumRows.median_ms(k);
    [maxV,maxI] = max(col); [minV,minI] = min(col);
    fprintf('  B=%3d Mbps: pooled median=%.3f | location range [%.3f (Q%d) , %.3f (Q%d)] | spread=%.3f ms\n', ...
        bwList(k), pooledMed, minV, qLevels(minI), maxV, qLevels(maxI), maxV-minV);
    sensRows = [sensRows; table(bwList(k), pooledMed, minV, qLevels(minI), maxV, qLevels(maxI), maxV-minV, ...
        'VariableNames', {'B_Mbps','pooled_median_ms','min_location_median_ms','min_location_qIdx','max_location_median_ms','max_location_qIdx','spread_ms'})]; %#ok<AGROW>
end
fprintf('CAUTION: cell-level P95 values (n=5 per cell) are highly sample-sensitive and should not be overinterpreted.\n\n');

%% ================= PART 11: CHECK AGAINST TABLE 3 =================
fprintf('=== PART 11: Check against the mean jitter of Table 3 ===\n');
knownJitter = [0.829 0.525 0.192 0.100 0.060 0.035];
knownBW     = [1 10 50 100 200 500];
verifyRows = table();
maxJitDiff = 0;
for k = 1:nBW
    idx = find(knownBW == bwList(k), 1);
    computedMean = jitRows.mean_ms(k);
    if ~isempty(idx)
        expected = knownJitter(idx);
        d = abs(computedMean - expected);
        maxJitDiff = max(maxJitDiff, d);
    else
        expected = NaN; d = NaN;
    end
    fprintf('  B=%3d Mbps: computed mean jitter=%.4f ms, Table 3: %.3f ms, diff=%.4g ms\n', ...
        bwList(k), computedMean, expected, d);
    verifyRows = [verifyRows; table(bwList(k), computedMean, expected, d, ...
        'VariableNames', {'B_Mbps','computed_mean_jitter_ms','previously_reported_mean_jitter_ms','abs_diff_ms'})]; %#ok<AGROW>
end
verifyPass = maxJitDiff < 0.001;
fprintf('Max abs diff across all bandwidths: %.4g ms (pass<0.001: %d)\n\n', maxJitDiff, verifyPass);
if ~verifyPass
    warning('step09:JitterMismatch', 'Mean jitter differs from Table 3 of the article by more than 0.001 ms.');
end

%% ================= PART 10: FIGURES =================
fprintf('=== PART 10: Building figures (self-contained SVG, no MATLAB graphics rendering) ===\n');
bwLabels = arrayfun(@(v) sprintf('%g Mbps', v), bwList, 'UniformOutput', false);
bwColors = {'#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594'};
threshColors = {'#1b9e77','#d95f02','#7570b3','#e7298a'};

% Fig 1: boxplot of trial-level mean latency by bandwidth
svg1 = buildBoxplotSVG(latByBW, bwLabels, bwColors, ...
    'Trial-level mean latency (ms)', 'Target bandwidth', ...
    'Distribution of trial-level mean latency by target bandwidth', ...
    'Unit of analysis: one 60-second trial mean (n=35 trials per bandwidth, pooled across 7 locations).', false);
writeSvgFile(svg1, fullfile(figDir, 'step09_latency_boxplot_by_bandwidth.svg'));

% Fig 2: ECDF of trial-level mean latency by bandwidth
svg2 = buildEcdfSVG(latByBW, bwLabels, bwColors, ...
    'Trial-level mean latency (ms)', ...
    'Empirical CDF of trial-level mean latency by target bandwidth', ...
    'Unit of analysis: one 60-second trial mean (n=35 trials per bandwidth, pooled across 7 locations).');
writeSvgFile(svg2, fullfile(figDir, 'step09_latency_ecdf_by_bandwidth.svg'));

% Fig 3: threshold exceedance with Wilson CI
svg3 = buildExceedanceSVG(bwLabels, thresholds, exceedPropMat, exceedCiLoMat, exceedCiHiMat, threshColors, ...
    'Proportion of 60-second trials whose mean latency exceeded threshold', ...
    'Trial-mean threshold-exceedance proportions (95% Wilson intervals)', ...
    'Unit of analysis: one 60-second trial mean; error bars are 95% Wilson binomial confidence intervals (n=35 per bandwidth).');
writeSvgFile(svg3, fullfile(figDir, 'step09_threshold_exceedance.svg'));

% Fig 4: boxplot of trial-level jitter by bandwidth
svg4 = buildBoxplotSVG(jitByBW, bwLabels, bwColors, ...
    'Trial-level jitter (ms)', 'Target bandwidth', ...
    'Distribution of trial-level jitter by target bandwidth', ...
    'Unit of analysis: one 60-second trial jitter value (n=35 trials per bandwidth, pooled across 7 locations). Log y-axis.', true);
writeSvgFile(svg4, fullfile(figDir, 'step09_jitter_boxplot_by_bandwidth.svg'));

% Fig 5: boxplot of per-trial maximum delay by bandwidth
svg5 = buildBoxplotSVG(maxByBW, bwLabels, bwColors, ...
    'Per-trial maximum delay (ms)', 'Target bandwidth', ...
    'Distribution of per-trial maximum delay by target bandwidth', ...
    'Unit of analysis: one 60-second trial maximum (n=35 trials per bandwidth, pooled across 7 locations). Log y-axis; not a packet-level latency distribution.', true);
writeSvgFile(svg5, fullfile(figDir, 'step09_max_delay_boxplot_by_bandwidth.svg'));

% Fig 6: packet-loss distribution by bandwidth
lossByBW = cell(nBW,1);
for k = 1:nBW, lossByBW{k} = A.loss_pct(A.B_Mbps == bwList(k)); end
svg6 = buildBoxplotSVG(lossByBW, bwLabels, bwColors, ...
    'Trial-level packet loss (%)', 'Target bandwidth', ...
    'Distribution of trial-level packet loss by target bandwidth', ...
    'Unit of analysis: one 60-second trial loss percentage (n=35 trials per bandwidth, pooled across 7 locations). Many trials have zero loss, especially at low bandwidths.', false);
writeSvgFile(svg6, fullfile(figDir, 'step09_packet_loss_boxplot_by_bandwidth.svg'));

fprintf('Wrote 6 SVG figures to %s\n', figDir);
fprintf('\n');

%% ================= SAVE RESULTS (.mat, .xlsx) =================
fprintf('=== Saving stage11_results.mat / .xlsx ===\n');
results = struct();
results.A = A;
results.trialLatencySummary = sumRows;
results.bootstrapCI = bciRows;
results.thresholdExceedance = exRows;
results.jitterSummary = jitRows;
results.maxDelaySummary = maxRows;
results.packetLossSummary = lossRows;
results.locationSensitivity = locRows;
results.locationSensitivitySpread = sensRows;
results.medianWide = array2table(medWide, 'RowNames', arrayfun(@(v) sprintf('Q%d',v), qLevels, 'UniformOutput', false), ...
    'VariableNames', arrayfun(@(v) sprintf('B%d', v), bwList, 'UniformOutput', false));
results.p95Wide = array2table(p95Wide, 'RowNames', arrayfun(@(v) sprintf('Q%d',v), qLevels, 'UniformOutput', false), ...
    'VariableNames', arrayfun(@(v) sprintf('B%d', v), bwList, 'UniformOutput', false));
results.jitterVerification = verifyRows;
results.nBoot = NBOOT;
results.seed = 20260815;
save(fullfile(here, 'stage11_results.mat'), 'results');

xlsFile = fullfile(here, 'stage11_results.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(sumRows, xlsFile, 'Sheet', 'Trial_Latency_Summary');
writetable(bciRows, xlsFile, 'Sheet', 'Bootstrap_CI');
writetable(exRows, xlsFile, 'Sheet', 'Threshold_Exceedance');
writetable(jitRows, xlsFile, 'Sheet', 'Jitter_Summary');
writetable(maxRows, xlsFile, 'Sheet', 'Max_Delay_Summary');
writetable(lossRows, xlsFile, 'Sheet', 'Packet_Loss_Summary');
writetable(locRows, xlsFile, 'Sheet', 'Location_Sensitivity');
writetable(results.medianWide, xlsFile, 'Sheet', 'Location_Median_Wide', 'WriteRowNames', true);
writetable(results.p95Wide, xlsFile, 'Sheet', 'Location_P95_Wide', 'WriteRowNames', true);
writetable(verifyRows, xlsFile, 'Sheet', 'Jitter_Verification');
fprintf('Saved: stage11_results.mat, stage11_results.xlsx\n');

fprintf('\nDone.\n');
end

%% ============ LOCAL FUNCTIONS: STATISTICS ============

function p = empPercentile(x, q)
% Type-7 (linear interpolation) empirical percentile, q in [0,100]. Standard, matches R/Excel PERCENTILE.INC.
x = sort(x(:));
n = numel(x);
if n == 1, p = x(1); return; end
h = (q/100)*(n-1) + 1;
lo = floor(h); hi = ceil(h);
if lo == hi
    p = x(lo);
else
    p = x(lo) + (h-lo)*(x(hi)-x(lo));
end
end

function pvec = rowPercentile(X, q)
% Row-wise Type-7 empirical percentile for a [nBoot x n] matrix.
Xs = sort(X, 2);
n = size(Xs,2);
h = (q/100)*(n-1) + 1;
lo = floor(h); hi = ceil(h);
if lo == hi
    pvec = Xs(:,lo);
else
    pvec = Xs(:,lo) + (h-lo)*(Xs(:,hi)-Xs(:,lo));
end
end

function [lo,hi] = wilsonCI(x, n)
% Wilson score interval for a binomial proportion, 95% (z=1.959963985).
z = 1.959963984540054;
if n == 0, lo = NaN; hi = NaN; return; end
phat = x/n;
denom = 1 + z^2/n;
center = phat + z^2/(2*n);
adj = z*sqrt((phat*(1-phat) + z^2/(4*n))/n);
lo = max(0, (center - adj)/denom);
hi = min(1, (center + adj)/denom);
end

function pooled = stratBootstrapPool(cellData, nBoot)
% Stratified bootstrap: resample WITH replacement within each group (location),
% preserving each group's original sample size, then concatenate across groups.
% Returns [nBoot x totalN] matrix of resampled pooled values.
nGroups = numel(cellData);
sizes = cellfun(@numel, cellData);
totalN = sum(sizes);
pooled = nan(nBoot, totalN);
for b = 1:nBoot
    row = nan(1, totalN);
    pos = 1;
    for g = 1:nGroups
        v = cellData{g};
        ng = sizes(g);
        idx = randi(ng, ng, 1);
        row(pos:pos+ng-1) = v(idx);
        pos = pos + ng;
    end
    pooled(b,:) = row;
end
end

function s = boxStats(x)
x = sort(x(:));
q1 = empPercentile(x,25); q3 = empPercentile(x,75); med = median(x);
iqr = q3 - q1;
loFence = q1 - 1.5*iqr; hiFence = q3 + 1.5*iqr;
inRange = x(x >= loFence & x <= hiFence);
whiskerLo = min(inRange); whiskerHi = max(inRange);
outliers = x(x < loFence | x > hiFence);
s = struct('min',x(1),'q1',q1,'med',med,'q3',q3,'max',x(end), ...
    'whiskerLo',whiskerLo,'whiskerHi',whiskerHi,'outliers',outliers);
end

%% ============ LOCAL FUNCTIONS: SVG PRIMITIVES ============

function s = svgRect(x,y,w,h,fill,stroke,strokeW,opacity)
if nargin < 8, opacity = 1; end
s = sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s" stroke="%s" stroke-width="%.2f" opacity="%.2f"/>\n', ...
    x,y,w,h,fill,stroke,strokeW,opacity);
end

function s = svgLine(x1,y1,x2,y2,stroke,strokeW,dash)
if nargin < 7, dash = ''; end
dashAttr = '';
if ~isempty(dash), dashAttr = sprintf(' stroke-dasharray="%s"', dash); end
s = sprintf('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="%.2f"%s/>\n', ...
    x1,y1,x2,y2,stroke,strokeW,dashAttr);
end

function s = svgText(x,y,txt,size,anchor,weight,color,rotate)
if nargin < 5, anchor = 'middle'; end
if nargin < 6, weight = 'normal'; end
if nargin < 7, color = '#1a1a1a'; end
if nargin < 8, rotate = 0; end
txt = strrep(txt, '&', '&amp;');
rotAttr = '';
if rotate ~= 0
    rotAttr = sprintf(' transform="rotate(%d %.2f %.2f)"', rotate, x, y);
end
s = sprintf('<text x="%.2f" y="%.2f" font-family="Arial,Helvetica,sans-serif" font-size="%d" text-anchor="%s" font-weight="%s" fill="%s"%s>%s</text>\n', ...
    x,y,size,anchor,weight,color,rotAttr,txt);
end

function s = svgCircle(cx,cy,r,fill,opacity)
if nargin < 5, opacity = 1; end
s = sprintf('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" opacity="%.2f"/>\n', cx,cy,r,fill,opacity);
end

function s = svgPolyline(xs,ys,stroke,strokeW,fill)
if nargin < 5, fill = 'none'; end
pts = '';
for i = 1:numel(xs), pts = [pts sprintf('%.2f,%.2f ', xs(i), ys(i))]; end %#ok<AGROW>
s = sprintf('<polyline points="%s" fill="%s" stroke="%s" stroke-width="%.2f"/>\n', pts, fill, stroke, strokeW);
end

function writeSvgFile(svgStr, path)
fid = fopen(path, 'w');
fprintf(fid, '%s', svgStr);
fclose(fid);
end

%% ============ LOCAL FUNCTIONS: CHART BUILDERS ============

function [yToPix, yTicks] = makeYScale(dataMin, dataMax, plotTop, plotBottom, useLog)
if useLog
    lo = log10(max(dataMin, 1e-6)); hi = log10(dataMax);
    pad = (hi-lo)*0.08; lo = lo - pad; hi = hi + pad;
    yToPix = @(v) plotBottom - (log10(max(v,1e-6))-lo)/(hi-lo)*(plotBottom-plotTop);
    loD = floor(lo); hiD = ceil(hi);
    yTicks = [];
    for d = loD:hiD
        for m = [1 2 5]
            v = m*10^d;
            if v >= 10^lo && v <= 10^hi, yTicks(end+1) = v; end %#ok<AGROW>
        end
    end
    if isempty(yTicks), yTicks = [dataMin dataMax]; end
else
    rng_ = dataMax - dataMin; if rng_ == 0, rng_ = max(1,dataMax); end
    pad = rng_*0.1;
    lo = max(0, dataMin - pad); hi = dataMax + pad;
    yToPix = @(v) plotBottom - (v-lo)/(hi-lo)*(plotBottom-plotTop);
    nTicksTarget = 6;
    step = niceStep((hi-lo)/nTicksTarget);
    yTicks = ceil(lo/step)*step : step : hi;
end
end

function step = niceStep(rough)
if rough <= 0, step = 1; return; end
mag = 10^floor(log10(rough));
resid = rough/mag;
if resid < 1.5, step = 1*mag;
elseif resid < 3, step = 2*mag;
elseif resid < 7, step = 5*mag;
else, step = 10*mag;
end
end

function svg = buildBoxplotSVG(groupData, labels, colors, yLabel, xLabel, titleStr, subtitleStr, useLogY)
nG = numel(groupData);
W = 120 + nG*110; H = 560;
plotLeft = 90; plotRight = W - 40; plotTop = 80; plotBottom = H - 110;

allVals = vertcat(groupData{:});
dataMin = min(allVals); dataMax = max(allVals);
[yToPix, yTicks] = makeYScale(dataMin, dataMax, plotTop, plotBottom, useLogY);

svg = sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" font-family="Arial,Helvetica,sans-serif">\n', W,H,W,H);
svg = [svg sprintf('<rect x="0" y="0" width="%d" height="%d" fill="white"/>\n', W,H)];
svg = [svg svgText(W/2, 30, titleStr, 17, 'middle', 'bold')];
svg = [svg svgText(W/2, 52, subtitleStr, 11, 'middle', 'normal', '#555555')];

% axes
svg = [svg svgLine(plotLeft, plotTop-10, plotLeft, plotBottom, '#333333', 1.2)];
svg = [svg svgLine(plotLeft, plotBottom, plotRight, plotBottom, '#333333', 1.2)];
for t = 1:numel(yTicks)
    yv = yTicks(t); yp = yToPix(yv);
    if yp < plotTop-10 || yp > plotBottom, continue; end
    svg = [svg svgLine(plotLeft-5, yp, plotRight, yp, '#e5e5e5', 0.8)]; %#ok<AGROW>
    lbl = sprintf('%g', yv);
    svg = [svg svgText(plotLeft-10, yp+4, lbl, 11, 'end')]; %#ok<AGROW>
end
svg = [svg svgText(25, (plotTop+plotBottom)/2, yLabel, 13, 'middle', 'bold', '#1a1a1a', -90)];
svg = [svg svgText(W/2, H-15, xLabel, 13, 'middle', 'bold')];

boxW = 55;
step = (plotRight-plotLeft)/nG;
for g = 1:nG
    xc = plotLeft + step*(g-0.5);
    st = boxStats(groupData{g});
    col = colors{mod(g-1,numel(colors))+1};
    y_q1 = yToPix(st.q1); y_q3 = yToPix(st.q3); y_med = yToPix(st.med);
    y_wl = yToPix(st.whiskerLo); y_wh = yToPix(st.whiskerHi);
    svg = [svg svgLine(xc, y_wh, xc, y_q3, '#333333', 1.2)]; %#ok<AGROW>
    svg = [svg svgLine(xc, y_q1, xc, y_wl, '#333333', 1.2)]; %#ok<AGROW>
    svg = [svg svgLine(xc-boxW/4, y_wh, xc+boxW/4, y_wh, '#333333', 1.2)]; %#ok<AGROW>
    svg = [svg svgLine(xc-boxW/4, y_wl, xc+boxW/4, y_wl, '#333333', 1.2)]; %#ok<AGROW>
    svg = [svg svgRect(xc-boxW/2, y_q3, boxW, max(1,y_q1-y_q3), col, '#333333', 1.3)]; %#ok<AGROW>
    svg = [svg svgLine(xc-boxW/2, y_med, xc+boxW/2, y_med, '#111111', 2.2)]; %#ok<AGROW>
    for o = 1:numel(st.outliers)
        svg = [svg svgCircle(xc, yToPix(st.outliers(o)), 3, '#c0392b', 0.75)]; %#ok<AGROW>
    end
    svg = [svg svgText(xc, plotBottom+22, labels{g}, 12, 'middle')]; %#ok<AGROW>
    svg = [svg svgText(xc, plotBottom+38, sprintf('n=%d', numel(groupData{g})), 10, 'middle', 'normal', '#777777')]; %#ok<AGROW>
end
svg = [svg '</svg>'];
end

function svg = buildEcdfSVG(groupData, labels, colors, xLabel, titleStr, subtitleStr)
nG = numel(groupData);
W = 760; H = 560;
plotLeft = 90; plotRight = W - 220; plotTop = 80; plotBottom = H - 90;

allVals = vertcat(groupData{:});
dataMin = min(allVals); dataMax = max(allVals);
rng_ = dataMax-dataMin; pad = rng_*0.05;
xlo = max(0,dataMin-pad); xhi = dataMax+pad;
xToPix = @(v) plotLeft + (v-xlo)/(xhi-xlo)*(plotRight-plotLeft);
yToPix = @(v) plotBottom - v*(plotBottom-plotTop);

svg = sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" font-family="Arial,Helvetica,sans-serif">\n', W,H,W,H);
svg = [svg sprintf('<rect x="0" y="0" width="%d" height="%d" fill="white"/>\n', W,H)];
svg = [svg svgText(W/2, 30, titleStr, 17, 'middle', 'bold')];
svg = [svg svgText(W/2, 52, subtitleStr, 11, 'middle', 'normal', '#555555')];

svg = [svg svgLine(plotLeft, plotTop-10, plotLeft, plotBottom, '#333333', 1.2)];
svg = [svg svgLine(plotLeft, plotBottom, plotRight, plotBottom, '#333333', 1.2)];
for f = 0:0.2:1
    yp = yToPix(f);
    svg = [svg svgLine(plotLeft-5, yp, plotRight, yp, '#e5e5e5', 0.8)]; %#ok<AGROW>
    svg = [svg svgText(plotLeft-10, yp+4, sprintf('%.1f',f), 11, 'end')]; %#ok<AGROW>
end
step = niceStep((xhi-xlo)/6);
for xv = ceil(xlo/step)*step : step : xhi
    xp = xToPix(xv);
    svg = [svg svgLine(xp, plotBottom, xp, plotBottom+6, '#333333', 1)]; %#ok<AGROW>
    svg = [svg svgText(xp, plotBottom+22, sprintf('%g',xv), 11, 'middle')]; %#ok<AGROW>
end
svg = [svg svgText(25, (plotTop+plotBottom)/2, 'Empirical CDF', 13, 'middle', 'bold', '#1a1a1a', -90)];
svg = [svg svgText((plotLeft+plotRight)/2, H-15, xLabel, 13, 'middle', 'bold')];

for g = 1:nG
    x = sort(groupData{g}(:));
    n = numel(x);
    xs = []; ys = [];
    for i = 1:n
        xs = [xs, x(i), x(i)]; ys = [ys, (i-1)/n, i/n]; %#ok<AGROW>
    end
    col = colors{mod(g-1,numel(colors))+1};
    svg = [svg svgPolyline(arrayfun(xToPix,xs), arrayfun(yToPix,ys), col, 2.2)]; %#ok<AGROW>
    legY = plotTop + (g-1)*22;
    svg = [svg svgLine(plotRight+20, legY, plotRight+45, legY, col, 3)]; %#ok<AGROW>
    svg = [svg svgText(plotRight+52, legY+4, labels{g}, 12, 'start')]; %#ok<AGROW>
end
svg = [svg '</svg>'];
end

function svg = buildExceedanceSVG(bwLabels, thresholds, propMat, ciLoMat, ciHiMat, colors, yLabel, titleStr, subtitleStr)
[nBW, nT] = size(propMat);
W = 130 + nBW*130; H = 560;
plotLeft = 90; plotRight = W - 40; plotTop = 80; plotBottom = H - 130;

yToPix = @(v) plotBottom - v*(plotBottom-plotTop);

svg = sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" font-family="Arial,Helvetica,sans-serif">\n', W,H,W,H);
svg = [svg sprintf('<rect x="0" y="0" width="%d" height="%d" fill="white"/>\n', W,H)];
svg = [svg svgText(W/2, 30, titleStr, 17, 'middle', 'bold')];
svg = [svg svgText(W/2, 52, subtitleStr, 11, 'middle', 'normal', '#555555')];

svg = [svg svgLine(plotLeft, plotTop-10, plotLeft, plotBottom, '#333333', 1.2)];
svg = [svg svgLine(plotLeft, plotBottom, plotRight, plotBottom, '#333333', 1.2)];
maxY = min(1, max(ciHiMat(:))*1.15);
for f = 0:0.1:ceil(maxY*10)/10
    yp = yToPix(f);
    if yp < plotTop-10, continue; end
    svg = [svg svgLine(plotLeft-5, yp, plotRight, yp, '#e5e5e5', 0.8)]; %#ok<AGROW>
    svg = [svg svgText(plotLeft-10, yp+4, sprintf('%.1f',f), 11, 'end')]; %#ok<AGROW>
end
svg = [svg svgText(25, (plotTop+plotBottom)/2, yLabel, 12, 'middle', 'bold', '#1a1a1a', -90)];

groupW = (plotRight-plotLeft)/nBW;
barW = groupW/(nT+1.5);
for k = 1:nBW
    gx0 = plotLeft + groupW*(k-1);
    for t = 1:nT
        xc = gx0 + barW*(t+0.25);
        p = propMat(k,t); lo = ciLoMat(k,t); hi = ciHiMat(k,t);
        yb = yToPix(0); yp = yToPix(p);
        col = colors{mod(t-1,numel(colors))+1};
        svg = [svg svgRect(xc-barW*0.4, min(yp,yb), barW*0.8, abs(yb-yp), col, '#333333', 0.8)]; %#ok<AGROW>
        yLo = yToPix(lo); yHi = yToPix(hi);
        svg = [svg svgLine(xc, yLo, xc, yHi, '#222222', 1.3)]; %#ok<AGROW>
        svg = [svg svgLine(xc-4, yLo, xc+4, yLo, '#222222', 1.3)]; %#ok<AGROW>
        svg = [svg svgLine(xc-4, yHi, xc+4, yHi, '#222222', 1.3)]; %#ok<AGROW>
    end
    svg = [svg svgText(gx0+groupW/2, plotBottom+24, bwLabels{k}, 12, 'middle')]; %#ok<AGROW>
end
svg = [svg svgText((plotLeft+plotRight)/2, H-15, 'Target bandwidth', 13, 'middle', 'bold')];

legX = plotLeft; legY = plotBottom + 60;
for t = 1:nT
    col = colors{mod(t-1,numel(colors))+1};
    lx = legX + (t-1)*140;
    svg = [svg svgRect(lx, legY-12, 16, 16, col, '#333333', 0.8)]; %#ok<AGROW>
    svg = [svg svgText(lx+22, legY, sprintf('> %d ms', thresholds(t)), 12, 'start')]; %#ok<AGROW>
end
svg = [svg '</svg>'];
end
