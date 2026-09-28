function step06_gamma_bandwidth_corr()
% STEP 06 - Location offsets versus location-specific throughput variability
% (Section 5.2).
%
% No radio-layer measurements (RSRP, RSRQ, SINR, MCS, retransmissions) are
% available, so the location-wise standard deviation of the achieved
% throughput is used as an indirect indicator of radio-channel variability
% and correlated with the location offsets gamma_q.
%
% gamma_q is taken from stage2_gamma_offsets.csv written by step02_fit_model
% (column gamma_Stage2: the full-data fit with the achieved throughput).
%
% Reads data/ and the output of step02_fit_model; writes outputs to the code folder.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);

%% ---- Step 1: gamma_q (full-sample, achieved-bandwidth fit) ----
gT = readtable(fullfile(here, 'stage2_gamma_offsets.csv'));
qIdx = gT.qIdx(:);
gamma = gT.gamma_Stage2(:);
fprintf('Loaded gamma_q (full-data fit, step 02):\n');
disp(table(qIdx, gamma));

%% ---- Step 2: per-location SD of achieved bandwidth ----
% (a) AS SPECIFIED: pooled across all trials AND all 6 target-bandwidth
%     levels (1/10/50/100/200/500 Mbps) into one SD per location.
% (b) DIAGNOSTIC CHECK: because the design is balanced (every location has
%     the SAME 6 target levels x 5 trials), the between-level spread
%     (1 Mbps to ~524 Mbps) swamps any location-specific radio variability
%     in (a). (b) instead
%     computes SD across the 5 trials WITHIN each target level, then
%     averages across the 6 levels per location - this isolates trial-to-
%     trial fluctuation AT a fixed target rate, which is what the stated
%     rationale (SINR/MCS-driven throughput jitter) actually describes.
T = readtable(fullfile(root, 'data', 'latency_with_achieved_bw_210.csv'));
nLoc = numel(qIdx);
levels = unique(T.B_Mbps);
bwSD_pooled = zeros(nLoc,1);
bwSD_withinLevel = zeros(nLoc,1);
nRows = zeros(nLoc,1);
for i = 1:nLoc
    mask = T.qIdx == qIdx(i);
    bwSD_pooled(i) = std(T.achieved_Mbps_server(mask));
    nRows(i) = sum(mask);
    lvlSD = zeros(numel(levels),1);
    for k = 1:numel(levels)
        m2 = mask & T.B_Mbps == levels(k);
        lvlSD(k) = std(T.achieved_Mbps_server(m2));
    end
    bwSD_withinLevel(i) = mean(lvlSD);
end
fprintf('\nPer-location achieved-bandwidth SD:\n');
disp(table(qIdx, bwSD_pooled, bwSD_withinLevel, nRows));

%% ---- Step 3: Pearson correlation (n=7), both metrics ----
[r_pooled, p_pooled] = corr(bwSD_pooled, gamma, 'type', 'Pearson');
[r_within, p_within] = corr(bwSD_withinLevel, gamma, 'type', 'Pearson');
fprintf('\n(a) Pooled SD:       Pearson r = %.4f, p = %.4f (n=%d)\n', r_pooled, p_pooled, nLoc);
fprintf('(b) Within-level SD: Pearson r = %.4f, p = %.4f (n=%d)\n', r_within, p_within, nLoc);

%% ---- Save outputs ----
resultTable = table(qIdx, gamma, bwSD_pooled, bwSD_withinLevel, nRows, ...
    'VariableNames', {'qIdx','gamma_q','bandwidth_SD_pooled_Mbps','bandwidth_SD_withinLevel_Mbps','n_rows'});
writetable(resultTable, fullfile(here, 'stage8_gamma_vs_bandwidthSD.csv'));

results = struct();
results.table = resultTable;
results.r_pooled = r_pooled;   results.p_pooled = p_pooled;
results.r_within = r_within;   results.p_within = p_within;
results.n = nLoc;
save(fullfile(here, 'stage8_results.mat'), 'results');

xlsFile = fullfile(here, 'stage8_results.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(resultTable, xlsFile, 'Sheet', 'GammaVsBandwidthSD');

fprintf('\nSaved: stage8_gamma_vs_bandwidthSD.csv, stage8_results.xlsx, stage8_results.mat\n');
fprintf('Done.\n');

end
