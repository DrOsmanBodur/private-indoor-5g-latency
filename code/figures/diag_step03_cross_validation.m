function diag_step03_cross_validation()
% Diagnostic figures for step 03 (out-of-sample validation).
%   1. MAE, RMSE and R^2 on the training trials (1-3) and the test trials (4-5)
%   2. leave-one-location-out R^2 and MAE per held-out location
%   3. predicted versus measured latency on the test trials
%   4. R^2 in-sample, on the test trials and pooled over the LOLO folds
%
% Inputs:  stage3_results.mat and stage2_results.mat written by steps 03 and 02,
%          data/latency_with_achieved_bw_210.csv
% Outputs: PNG files in figures/diagnostics/

here = fileparts(mfilename('fullpath'));
codeDir = fileparts(here);
root = fileparts(codeDir);
figDir = fullfile(here, 'diagnostics');
if ~exist(figDir, 'dir'), mkdir(figDir); end

S3 = load(fullfile(codeDir, 'stage3_results.mat')); s3 = S3.results;
S2 = load(fullfile(codeDir, 'stage2_results.mat')); s2 = S2.results;

T = readtable(fullfile(root, 'data', 'latency_with_achieved_bw_210.csv'));
qIdx  = T.qIdx(:); trial = T.trial(:);
b     = T.achieved_Mbps_server(:); L = T.latency_ms(:);
testMask = trial >= 4;

%% ============= FIGURE 1: train vs test metrics =============
f = figure('Visible','off','Position',[100 100 950 380]);
tl = tiledlayout(f, 1, 3, 'Padding','compact','TileSpacing','compact');

nexttile(tl);
bar([s3.trainMetrics.MAE, s3.testMetrics.MAE]);
set(gca,'XTickLabel',{'Train (in-sample)','Test (out-of-sample)'});
title('MAE (ms)'); grid on; box on;

nexttile(tl);
bar([s3.trainMetrics.RMSE, s3.testMetrics.RMSE]);
set(gca,'XTickLabel',{'Train (in-sample)','Test (out-of-sample)'});
title('RMSE (ms)'); grid on; box on;

nexttile(tl);
bar([s3.trainMetrics.R2, s3.testMetrics.R2]);
set(gca,'XTickLabel',{'Train (in-sample)','Test (out-of-sample)'});
ylim([0.7 0.9]);
title('R^2'); grid on; box on;

title(tl, 'Trial-based split: training trials 1-3, test trials 4-5');
saveas(f, fullfile(figDir, 'step03_train_test_metrics.png'));
close(f);

%% ============= FIGURE 2: LOLO per-location R2 and MAE =============
lolo = s3.loloPerFold;
f = figure('Visible','off','Position',[100 100 950 550]);
tl = tiledlayout(f, 1, 2, 'Padding','compact','TileSpacing','compact');

nexttile(tl);
bb = bar(lolo.heldOut_qIdx, lolo.R2);
bb.FaceColor = [0.16 0.47 0.84];
yline(s2.R2, '--', sprintf('in-sample R^2 = %.3f', s2.R2), 'Color',[0.92 0.41 0.20], 'LineWidth',1.5);
yline(s3.loloPooled.R2, ':', sprintf('LOLO pooled R^2 = %.3f', s3.loloPooled.R2), 'Color',[0.11 0.69 0.48], 'LineWidth',1.5);
set(gca,'XTick',1:7,'XTickLabel',arrayfun(@(k)sprintf('Q%d',k),1:7,'UniformOutput',false));
title('R^2 per held-out location'); ylim([0 1]); grid on; box on;

nexttile(tl);
bb2 = bar(lolo.heldOut_qIdx, lolo.MAE);
bb2.FaceColor = [0.16 0.47 0.84];
yline(s2.MAE, '--', sprintf('in-sample MAE = %.2f', s2.MAE), 'Color',[0.92 0.41 0.20], 'LineWidth',1.5);
set(gca,'XTick',1:7,'XTickLabel',arrayfun(@(k)sprintf('Q%d',k),1:7,'UniformOutput',false));
title('MAE (ms) per held-out location'); grid on; box on;

title(tl, 'Leave-one-location-out validation');
saveas(f, fullfile(figDir, 'step03_lolo_per_location.png'));
close(f);

%% ============= FIGURE 3: predicted vs measured, test trials =============
pA = s3.trainParams(1:4);
gA = s3.trainGamma;
qTest = round(qIdx(testMask)); qTest(qTest<1)=1; qTest(qTest>7)=7;
predTest = pA(1) + pA(2).*b(testMask) + pA(3)./(pA(4) - b(testMask)) + gA(qTest)';
actualTest = L(testMask);

f = figure('Visible','off','Position',[100 100 700 700]);
scatter(actualTest, predTest, 36, [0.16 0.47 0.84], 'filled', 'MarkerFaceAlpha',0.7);
hold on;
lims = [0, max([actualTest; predTest])*1.05];
plot(lims, lims, 'k--', 'LineWidth', 1.3);
xlim(lims); ylim(lims); axis square;
xlabel('Measured latency (ms), test trials 4-5');
ylabel('Predicted latency (ms)');
title({'Predicted versus measured latency, test trials (n = 84)', sprintf('MAE = %.2f ms, RMSE = %.2f ms, R^2 = %.3f', s3.testMetrics.MAE, s3.testMetrics.RMSE, s3.testMetrics.R2)});
grid on; box on;
saveas(f, fullfile(figDir, 'step03_predicted_vs_measured_test.png'));
close(f);

%% ============= FIGURE 4: R2 summary =============
f = figure('Visible','off','Position',[100 100 800 550]);
labels = {'In-sample (n = 210)','Test trials (n = 84)','LOLO pooled (n = 210)'};
vals = [s2.R2, s3.testMetrics.R2, s3.loloPooled.R2];
bb3 = bar(vals);
bb3.FaceColor = 'flat';
bb3.CData = [0.92 0.41 0.20; 0.16 0.47 0.84; 0.11 0.69 0.48];
set(gca, 'XTickLabel', labels);
ylim([0.75 0.86]);
ylabel('R^2');
title('R^2 in-sample and in the two out-of-sample validations');
grid on; box on;
for i=1:3
    text(i, vals(i)+0.003, sprintf('%.3f', vals(i)), 'HorizontalAlignment','center', 'FontWeight','bold');
end
saveas(f, fullfile(figDir, 'step03_r2_summary.png'));
close(f);

fprintf('Figures written to %s\n', figDir);
end
