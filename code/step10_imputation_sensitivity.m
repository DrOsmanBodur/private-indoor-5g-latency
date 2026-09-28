function step10_imputation_sensitivity()
% STEP 10 - Sensitivity of the 50 Mbps trial-level latency statistics to the
% three replaced runs (Q1, 50 Mbps, trials 3-5; Section 4.3.1).
%
% The one-way delay of these runs was invalid (32-bit counter wrap-around)
% and was replaced by the mean of trials 1-2 of the same cell. This step
% concerns only the trial-level MEAN LATENCY statistics of step 09; jitter,
% maximum delay and packet loss of these runs are valid and are not
% re-examined here.
%
% Method: mean, median, P90, P95 and the 10/20/30/50 ms threshold exceedance
% (with Wilson 95% CIs) at 50 Mbps are recomputed for
%   n=35 (full pooled sample, as used in step 09)
%   n=32 (the same sample without the 3 replaced runs)
% and the differences are reported.
%
% Writes outputs to the code folder (no data files modified).

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
rng(20260815, 'twister');
NBOOT = 10000;

mainFile = fullfile(root, 'data', 'latency_210.csv');
suppFile = fullfile(root, 'data', 'stage5_trial_level.csv');
Tmain = readtable(mainFile);
Tsupp = readtable(suppFile);
assert(isequal(Tmain.qIdx, Tsupp.qIdx) && isequal(Tmain.B_Mbps, Tsupp.B_Mbps) && isequal(Tmain.trial, Tsupp.trial), ...
    'Row order mismatch between latency_210.csv and stage5_trial_level.csv.');

B_TARGET = 50;
mask50 = Tmain.B_Mbps == B_TARGET;
qIdx50 = Tmain.qIdx(mask50);
trial50 = Tmain.trial(mask50);
lat50   = Tmain.latency_ms(mask50);
corrupt50 = Tsupp.was_corrupted_and_fixed(mask50) == 1;

assert(numel(lat50) == 35, 'Expected 35 trials at B=50 Mbps, got %d', numel(lat50));
assert(sum(corrupt50) == 3, 'Expected exactly 3 mean-filled rows at B=50 Mbps, got %d', sum(corrupt50));

fprintf('=== 50 Mbps: sensitivity to the three replaced runs ===\n');
fprintf('Mean-filled rows being tested (qIdx, trial, latency_ms):\n');
disp(table(qIdx50(corrupt50), trial50(corrupt50), lat50(corrupt50), 'VariableNames', {'qIdx','trial','latency_ms'}));

full35 = lat50;                 % n=35, includes the 3 mean-filled rows
clean32 = lat50(~corrupt50);    % n=32, mean-filled rows excluded
qIdxFull35 = qIdx50;
qIdxClean32 = qIdx50(~corrupt50);

fprintf('n(full)=%d, n(excl. mean-filled)=%d\n\n', numel(full35), numel(clean32));

%% ---- descriptive point estimates ----
descRows = table();
[descRows] = addDescRow(descRows, 'n=35 (full, incl. 3 mean-filled)', full35);
[descRows] = addDescRow(descRows, 'n=32 (excl. 3 mean-filled)', clean32);
disp(descRows);

%% ---- bootstrap CIs (stratified by location within each sample) ----
% n=35 sample: 7 locations, Q1 has 5 trials (3 mean-filled), others have 5 genuine.
% n=32 sample: 7 locations, Q1 has only 2 trials (unbalanced - documented explicitly).
qLevels = unique(qIdx50);
cellData35 = cell(numel(qLevels),1);
cellData32 = cell(numel(qLevels),1);
for i = 1:numel(qLevels)
    cellData35{i} = full35(qIdxFull35 == qLevels(i));
    cellData32{i} = clean32(qIdxClean32 == qLevels(i));
end
fprintf('n=32 sample location cell sizes: %s (Q1 reduced from 5 to 2 - documented, not silently rebalanced)\n\n', ...
    mat2str(cellfun(@numel, cellData32)'));

pooled35 = stratBootstrapPool(cellData35, NBOOT);
pooled32 = stratBootstrapPool(cellData32, NBOOT);

ciRows = table();
stats = {'mean','median','P90','P95'};
for s = 1:numel(stats)
    [pt35, lo35, hi35] = bootStat(full35, pooled35, stats{s});
    [pt32, lo32, hi32] = bootStat(clean32, pooled32, stats{s});
    row = table(string(stats{s}), pt35, lo35, hi35, pt32, lo32, hi32, pt32-pt35, ...
        'VariableNames', {'statistic','n35_point_ms','n35_CI95_lower_ms','n35_CI95_upper_ms', ...
                           'n32_point_ms','n32_CI95_lower_ms','n32_CI95_upper_ms','diff_n32_minus_n35_ms'});
    ciRows = [ciRows; row]; %#ok<AGROW>
end
disp(ciRows);

%% ---- threshold exceedance comparison ----
thresholds = [10 20 30 50];
exRows = table();
for t = 1:numel(thresholds)
    thr = thresholds(t);
    n35 = numel(full35); c35 = sum(full35 > thr); p35 = c35/n35; [lo35,hi35] = wilsonCI(c35,n35);
    n32 = numel(clean32); c32 = sum(clean32 > thr); p32 = c32/n32; [lo32,hi32] = wilsonCI(c32,n32);
    row = table(thr, n35, c35, p35, lo35, hi35, n32, c32, p32, lo32, hi32, p32-p35, ...
        'VariableNames', {'threshold_ms','n35','n_exceed_35','prop_35','CI35_lower','CI35_upper', ...
                           'n32','n_exceed_32','prop_32','CI32_lower','CI32_upper','diff_prop_32_minus_35'});
    exRows = [exRows; row]; %#ok<AGROW>
end
disp(exRows);

%% ---- verdict ----
maxPointDiff = max(abs(ciRows.n32_point_ms - ciRows.n35_point_ms));
maxPropDiff = max(abs(exRows.diff_prop_32_minus_35));
fprintf('\nMax abs difference in point estimates (mean/median/P90/P95), n=32 vs n=35: %.4f ms\n', maxPointDiff);
fprintf('Max abs difference in exceedance proportions, n=32 vs n=35: %.4f (%.1f percentage points)\n', maxPropDiff, maxPropDiff*100);

stable = (maxPointDiff < 0.5) && (maxPropDiff < 0.10);
fprintf('Conclusion stability (point-estimate diff < 0.5 ms AND exceedance-proportion diff < 10pp): %d\n\n', stable);

%% ---- save outputs ----
writetable(descRows, fullfile(here, 'stage11_addendum_b50_descriptive_comparison.csv'));
writetable(ciRows, fullfile(here, 'stage11_addendum_b50_bootstrap_comparison.csv'));
writetable(exRows, fullfile(here, 'stage11_addendum_b50_exceedance_comparison.csv'));

xlsFile = fullfile(here, 'stage11_addendum_b50_sensitivity.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(descRows, xlsFile, 'Sheet', 'Descriptive_Comparison');
writetable(ciRows, xlsFile, 'Sheet', 'Bootstrap_Comparison');
writetable(exRows, xlsFile, 'Sheet', 'Exceedance_Comparison');

results = struct('descRows', descRows, 'ciRows', ciRows, 'exRows', exRows, ...
    'full35', full35, 'clean32', clean32, 'maxPointDiff', maxPointDiff, 'maxPropDiff', maxPropDiff, 'stable', stable);
save(fullfile(here, 'stage11_addendum_b50_sensitivity.mat'), 'results');

fprintf('Saved: stage11_addendum_b50_descriptive_comparison.csv, stage11_addendum_b50_bootstrap_comparison.csv,\n');
fprintf('       stage11_addendum_b50_exceedance_comparison.csv, stage11_addendum_b50_sensitivity.xlsx/.mat\n');
fprintf('Done.\n');
end

%% ============ LOCAL FUNCTIONS ============

function T = addDescRow(T, label, x)
q1 = empPercentile(x,25); q3 = empPercentile(x,75);
row = table(string(label), numel(x), mean(x), std(x), median(x), q1, q3, q3-q1, ...
    empPercentile(x,90), empPercentile(x,95), min(x), max(x), ...
    'VariableNames', {'sample','n','mean_ms','sd_ms','median_ms','Q1_ms','Q3_ms','IQR_ms','P90_ms','P95_ms','min_ms','max_ms'});
T = [T; row];
end

function [pt, lo, hi] = bootStat(x, pooled, statName)
switch statName
    case 'mean',   pt = mean(x);   bvals = mean(pooled,2);
    case 'median', pt = median(x); bvals = median(pooled,2);
    case 'P90',    pt = empPercentile(x,90); bvals = rowPercentile(pooled,90);
    case 'P95',    pt = empPercentile(x,95); bvals = rowPercentile(pooled,95);
end
lo = empPercentile(bvals,2.5); hi = empPercentile(bvals,97.5);
end

function p = empPercentile(x, q)
x = sort(x(:));
n = numel(x);
if n == 1, p = x(1); return; end
h = (q/100)*(n-1) + 1;
lo = floor(h); hi = ceil(h);
if lo == hi, p = x(lo);
else, p = x(lo) + (h-lo)*(x(hi)-x(lo));
end
end

function pvec = rowPercentile(X, q)
Xs = sort(X, 2);
n = size(Xs,2);
h = (q/100)*(n-1) + 1;
lo = floor(h); hi = ceil(h);
if lo == hi, pvec = Xs(:,lo);
else, pvec = Xs(:,lo) + (h-lo)*(Xs(:,hi)-Xs(:,lo));
end
end

function [lo,hi] = wilsonCI(x, n)
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
        if ng == 0, continue; end
        idx = randi(ng, ng, 1);
        row(pos:pos+ng-1) = v(idx);
        pos = pos + ng;
    end
    pooled(b,:) = row;
end
end
