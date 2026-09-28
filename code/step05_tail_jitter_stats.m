function step05_tail_jitter_stats()
% STEP 05 - Trial-level jitter, packet-loss and delay statistics.
%
% Uses the Server Report fields parsed in step 01
% (data/achieved_bandwidth_210.csv); the raw logs are not re-parsed.
% Per (location qIdx, target rate B_Mbps) cell and overall:
%   - jitter_ms    : iPerf UDP jitter, as reported in the Server Report
%   - loss_pct     : packet loss, as reported in the Server Report
%   - P95/P99 (Gaussian): latency_avg + z*latency_stdev of the same trial
%                    (z95 = 1.645, z99 = 2.326). This approximation assumes
%                    approximately normal per-datagram delay within a trial;
%                    per-datagram delays are not available.
%   - observed max : the maximum delay reported for the trial, shown next to
%                    the Gaussian P99 approximation
%
% Known data issue: the three runs Q1, 50 Mbps, trials 3-5 have one-way delay
% avg/min values invalidated by the 32-bit counter wrap-around. For these
% rows latency_avg/min are replaced by the mean of trials 1-2 of the same
% cell (as in data/latency_210.csv); their jitter, loss, max and stdev are
% valid and used as reported.
%
% Reads data/achieved_bandwidth_210.csv; writes outputs to the code folder.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);

T = readtable(fullfile(root, 'data', 'achieved_bandwidth_210.csv'));
n = height(T);
fprintf('Loaded %d rows (expected 210) from data/achieved_bandwidth_210.csv\n', n);

qIdx = T.qIdx; B = T.B_target_Mbps; trial = T.trial;
lat_avg = T.latency_avg_ms_raw; lat_min = T.latency_min_ms_raw;
lat_max = T.latency_max_ms_raw; lat_std = T.latency_stdev_ms_raw;
jitter  = T.jitter_ms_raw; loss = T.loss_pct_raw;
corrupt = logical(T.latency_corrupted_flag);

fprintf('Corrupted rows (overflow bug): %d\n', sum(corrupt));
disp(T(corrupt, {'qIdx','B_target_Mbps','trial','latency_avg_ms_raw','latency_min_ms_raw','source_file'}));

%% ============= FIX CORRUPTED ROWS: mean-fill avg/min from trials 1-2 of same cell =============
fixIdx = find(corrupt);
for k = 1:numel(fixIdx)
    i = fixIdx(k);
    sameCell = qIdx == qIdx(i) & B == B(i) & trial <= 2 & ~corrupt;
    if any(sameCell)
        lat_avg(i) = mean(lat_avg(sameCell));
        lat_min(i) = mean(lat_min(sameCell));
        fprintf('Row %d (qIdx=%d,B=%d,trial=%d): avg/min mean-filled from trials 1-2 -> avg=%.3f min=%.3f\n', ...
            i, qIdx(i), B(i), trial(i), lat_avg(i), lat_min(i));
    else
        warning('Row %d: no clean same-cell rows to mean-fill from.', i);
    end
end

%% ============= PER-TRIAL TAIL ESTIMATES =============
z95 = 1.6449; z99 = 2.3263;
p95_est = lat_avg + z95 .* lat_std;
p99_est = lat_avg + z99 .* lat_std;

trialTable = table(qIdx, B, trial, lat_avg, lat_std, lat_min, lat_max, ...
    jitter, loss, p95_est, p99_est, corrupt, ...
    'VariableNames', {'qIdx','B_Mbps','trial','latency_avg_ms','latency_stdev_ms', ...
    'latency_min_ms','latency_max_ms','jitter_ms','loss_pct','P95_est_ms','P99_est_ms','was_corrupted_and_fixed'});
writetable(trialTable, fullfile(here, 'stage5_trial_level.csv'));

%% ============= AGGREGATE BY (location, bandwidth) CELL =============
cells = unique([qIdx, B], 'rows');
nCells = size(cells,1);
c_qIdx = cells(:,1); c_B = cells(:,2);
c_jitter_mean = nan(nCells,1); c_jitter_sd = nan(nCells,1);
c_loss_mean = nan(nCells,1); c_loss_max = nan(nCells,1);
c_p95_mean = nan(nCells,1); c_p99_mean = nan(nCells,1);
c_obsmax_mean = nan(nCells,1); c_obsmax_max = nan(nCells,1);
c_lat_mean = nan(nCells,1);

for c = 1:nCells
    idx = qIdx==c_qIdx(c) & B==c_B(c);
    c_jitter_mean(c) = mean(jitter(idx)); c_jitter_sd(c) = std(jitter(idx));
    c_loss_mean(c) = mean(loss(idx)); c_loss_max(c) = max(loss(idx));
    c_p95_mean(c) = mean(p95_est(idx)); c_p99_mean(c) = mean(p99_est(idx));
    c_obsmax_mean(c) = mean(lat_max(idx)); c_obsmax_max(c) = max(lat_max(idx));
    c_lat_mean(c) = mean(lat_avg(idx));
end
cellTable = table(c_qIdx, c_B, c_lat_mean, c_jitter_mean, c_jitter_sd, ...
    c_loss_mean, c_loss_max, c_p95_mean, c_p99_mean, c_obsmax_mean, c_obsmax_max, ...
    'VariableNames', {'qIdx','B_Mbps','latency_mean_ms','jitter_mean_ms','jitter_sd_ms', ...
    'loss_pct_mean','loss_pct_max','P95_est_mean_ms','P99_est_mean_ms', ...
    'observed_max_mean_ms','observed_max_worst_ms'});
cellTable = sortrows(cellTable, {'qIdx','B_Mbps'});
writetable(cellTable, fullfile(here, 'stage5_cell_summary.csv'));

%% ============= OVERALL SUMMARY =============
fprintf('\n=== OVERALL (n=%d trials) ===\n', n);
fprintf('Jitter (ms):    mean=%.4f  median=%.4f  min=%.4f  max=%.4f\n', ...
    mean(jitter), median(jitter), min(jitter), max(jitter));
fprintf('Loss (%%):       mean=%.4f  median=%.4f  max=%.4f  (%d/%d trials with loss>0)\n', ...
    mean(loss), median(loss), max(loss), sum(loss>0), n);
fprintf('P95 est. (ms):  mean=%.4f  min=%.4f  max=%.4f\n', mean(p95_est), min(p95_est), max(p95_est));
fprintf('P99 est. (ms):  mean=%.4f  min=%.4f  max=%.4f\n', mean(p99_est), min(p99_est), max(p99_est));
fprintf('Observed max (ms): mean=%.4f  worst single trial=%.4f\n', mean(lat_max), max(lat_max));
fprintf('Mean gap (observed_max - P99_est): %.4f ms  (heavier-than-Gaussian tail, as expected)\n', ...
    mean(lat_max - p99_est));

%% ============= JITTER / LOSS vs BANDWIDTH (does variability grow near capacity?) =============
uB = unique(B);
fprintf('\n--- Jitter and loss by target bandwidth (pooled over all 7 locations) ---\n');
fprintf('%8s %12s %12s %12s\n', 'B(Mbps)', 'jitter_mean', 'loss_mean%%', 'P99est_mean');
byBw = table();
for i = 1:numel(uB)
    idx = B==uB(i);
    fprintf('%8d %12.4f %12.4f %12.4f\n', uB(i), mean(jitter(idx)), mean(loss(idx)), mean(p99_est(idx)));
    byBw = [byBw; table(uB(i), mean(jitter(idx)), std(jitter(idx)), mean(loss(idx)), mean(p99_est(idx)), ...
        'VariableNames', {'B_Mbps','jitter_mean_ms','jitter_sd_ms','loss_pct_mean','P99_est_mean_ms'})]; %#ok<AGROW>
end
writetable(byBw, fullfile(here, 'stage5_by_bandwidth.csv'));

fprintf('\n--- Jitter and loss by location (pooled over all 6 bandwidth levels) ---\n');
fprintf('%8s %12s %12s %12s\n', 'Q', 'jitter_mean', 'loss_mean%%', 'P99est_mean');
byLoc = table();
for q = 1:7
    idx = qIdx==q;
    fprintf('%8d %12.4f %12.4f %12.4f\n', q, mean(jitter(idx)), mean(loss(idx)), mean(p99_est(idx)));
    byLoc = [byLoc; table(q, mean(jitter(idx)), std(jitter(idx)), mean(loss(idx)), mean(p99_est(idx)), ...
        'VariableNames', {'qIdx','jitter_mean_ms','jitter_sd_ms','loss_pct_mean','P99_est_mean_ms'})]; %#ok<AGROW>
end
writetable(byLoc, fullfile(here, 'stage5_by_location.csv'));

%% ============= CORRELATION: jitter/loss vs proximity to capacity (b/C) =============
% Uses the achieved throughput and the fitted effective capacity
% C = 611.39 Mbps (Table 5) as the capacity
% reference, to connect this stage back to the queueing story: does
% variability (jitter) and loss grow as b -> C?
C_ref = 611.39;
b_ach = T.achieved_Mbps_server;
proximity = b_ach ./ C_ref;
rJitter = corr(proximity, jitter, 'rows','complete');
rLoss = corr(proximity, loss, 'rows','complete');
fprintf('\n--- Correlation with proximity to capacity (b_achieved / C=%.2f) ---\n', C_ref);
fprintf('corr(proximity, jitter)  = %.4f\n', rJitter);
fprintf('corr(proximity, loss_pct) = %.4f\n', rLoss);

%% ============= SAVE =============
results = struct();
results.n = n; results.trialTable = trialTable; results.cellTable = cellTable;
results.byBw = byBw; results.byLoc = byLoc;
results.overall.jitter_mean = mean(jitter); results.overall.jitter_median = median(jitter);
results.overall.loss_mean = mean(loss); results.overall.loss_max = max(loss);
results.overall.p95_mean = mean(p95_est); results.overall.p99_mean = mean(p99_est);
results.overall.observed_max_mean = mean(lat_max); results.overall.observed_max_worst = max(lat_max);
results.corr_proximity_jitter = rJitter; results.corr_proximity_loss = rLoss;
results.C_ref = C_ref;
save(fullfile(here, 'stage5_results.mat'), 'results');

xlsFile = fullfile(here, 'stage5_results.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(trialTable, xlsFile, 'Sheet', 'Trial level');
writetable(cellTable, xlsFile, 'Sheet', 'Cell summary');
writetable(byBw, xlsFile, 'Sheet', 'By bandwidth');
writetable(byLoc, xlsFile, 'Sheet', 'By location');

fprintf('\nSaved: stage5_trial_level.csv, stage5_cell_summary.csv, stage5_by_bandwidth.csv,\n');
fprintf('       stage5_by_location.csv, stage5_results.mat, stage5_results.xlsx\n');
fprintf('Done.\n');

end
