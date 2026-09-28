function diag_step05_jitter_loss()
% Diagnostic figures for step 05 (trial-level jitter, loss and delay).
%   1. mean jitter (+/- SD) per target rate
%   2. mean packet loss per target rate
%   3. Gaussian P99 approximation versus the observed maximum delay
%   4. jitter and packet loss versus proximity to the effective capacity b/C
%
% Inputs:  stage5_results.mat written by step05_tail_jitter_stats,
%          data/achieved_bandwidth_210.csv, data/stage5_trial_level.csv
% Outputs: PNG files in figures/diagnostics/

here = fileparts(mfilename('fullpath'));
codeDir = fileparts(here);
root = fileparts(codeDir);
figDir = fullfile(here, 'diagnostics');
if ~exist(figDir, 'dir'), mkdir(figDir); end

S5 = load(fullfile(codeDir, 'stage5_results.mat')); s5 = S5.results;
byBw = s5.byBw;       % B_Mbps, jitter_mean_ms, jitter_sd_ms, loss_pct_mean, P99_est_mean_ms
cellT = s5.cellTable;  % qIdx, B_Mbps, ..., observed_max_mean_ms, observed_max_worst_ms

bwLevels = byBw.B_Mbps;
nB = numel(bwLevels);
bwLabels = arrayfun(@(k)sprintf('%dM',k),bwLevels,'UniformOutput',false);

%% ============= FIGURE 1: jitter vs target rate (mean +/- sd) =============
f = figure('Visible','off','Position',[100 100 900 550]);
errorbar(1:nB, byBw.jitter_mean_ms, byBw.jitter_sd_ms, 'o-', ...
    'Color',[0.16 0.47 0.84], 'MarkerFaceColor',[0.16 0.47 0.84], 'LineWidth',1.8, 'CapSize',8);
set(gca,'XTick',1:nB,'XTickLabel',bwLabels);
ylabel('Jitter (ms)'); xlabel('Target UDP rate');
title('Mean jitter (\pm SD) per target rate, pooled over the 7 locations');
grid on; box on;
saveas(f, fullfile(figDir, 'step05_jitter_vs_rate.png'));
close(f);

%% ============= FIGURE 2: packet loss vs target rate =============
f = figure('Visible','off','Position',[100 100 900 550]);
bar(1:nB, byBw.loss_pct_mean, 'FaceColor',[0.92 0.41 0.20], 'EdgeColor','none');
set(gca,'XTick',1:nB,'XTickLabel',bwLabels);
ylabel('Packet loss (%)'); xlabel('Target UDP rate');
title({'Mean packet loss per target rate', ...
       sprintf('correlation with proximity to capacity b/C: r = %.4f', s5.corr_proximity_loss)});
grid on; box on;
saveas(f, fullfile(figDir, 'step05_packet_loss_vs_rate.png'));
close(f);

%% ============= FIGURE 3: Gaussian P99 approximation vs observed maximum =============
obsMaxByBw = arrayfun(@(b) mean(cellT.observed_max_mean_ms(cellT.B_Mbps==b)), bwLevels);
ratio = obsMaxByBw ./ byBw.P99_est_mean_ms;

f = figure('Visible','off','Position',[100 100 950 600]);
x = 1:nB;
b1 = bar(x, [byBw.P99_est_mean_ms, obsMaxByBw]);
b1(1).FaceColor = [0.11 0.69 0.48];
b1(2).FaceColor = [0.75 0.10 0.10];
set(gca, 'YScale','log', 'XTick',1:nB,'XTickLabel',bwLabels);
ylabel('Delay (ms, log scale)'); xlabel('Target UDP rate');
legend({'Gaussian P99 approximation (avg + z\cdotstdev)','Observed maximum'}, 'Location','northwest');
title({'Gaussian P99 approximation versus observed maximum delay', ...
       sprintf('observed maximum / P99 approximation: %.0f to %.0f', min(ratio), max(ratio))});
grid on; box on;
saveas(f, fullfile(figDir, 'step05_p99_vs_observed_max.png'));
close(f);

%% ============= FIGURE 4: jitter/loss vs proximity to capacity =============
A = readtable(fullfile(root, 'data', 'achieved_bandwidth_210.csv'));
T5 = readtable(fullfile(root, 'data', 'stage5_trial_level.csv'));
proximity = A.achieved_Mbps_server ./ s5.C_ref;

f = figure('Visible','off','Position',[100 100 950 700]);
tl = tiledlayout(f, 2, 1, 'Padding','compact', 'TileSpacing','compact');

nexttile(tl);
scatter(proximity, T5.jitter_ms, 14, [0.16 0.47 0.84], 'filled', 'MarkerFaceAlpha', 0.45);
xlabel('Proximity to capacity (b_{achieved}/C)'); ylabel('Jitter (ms)');
title(sprintf('Jitter versus proximity to capacity (r = %.4f)', s5.corr_proximity_jitter));
grid on; box on;

nexttile(tl);
scatter(proximity, T5.loss_pct, 14, [0.92 0.41 0.20], 'filled', 'MarkerFaceAlpha', 0.45);
xlabel('Proximity to capacity (b_{achieved}/C)'); ylabel('Packet loss (%)');
title(sprintf('Packet loss versus proximity to capacity (r = %.4f)', s5.corr_proximity_loss));
grid on; box on;

title(tl, 'Jitter and packet loss versus proximity to capacity (n = 210 trials)');
saveas(f, fullfile(figDir, 'step05_proximity_scatter.png'));
close(f);

fprintf('Figures written to %s\n', figDir);
end
