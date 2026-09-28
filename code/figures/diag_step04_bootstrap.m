function diag_step04_bootstrap()
% Diagnostic figures for step 04 (parameter confidence intervals).
%   1. bootstrap distributions of L0, beta, theta and C with the point
%      estimate and the 95% percentile interval
%   2. location offsets with their bootstrap 95% intervals
%   3. width of the bootstrap and asymptotic 95% intervals of theta and C
%   4. bootstrap estimates of theta versus C (the theta-C ridge, Section 5.1)
%
% Input:   stage4_results.mat written by step04_bootstrap_ci
% Outputs: PNG files in figures/diagnostics/

here = fileparts(mfilename('fullpath'));
codeDir = fileparts(here);
figDir = fullfile(here, 'diagnostics');
if ~exist(figDir, 'dir'), mkdir(figDir); end

S4 = load(fullfile(codeDir, 'stage4_results.mat')); s4 = S4.results;
boot = s4.bootParams; % 1000x10: L0 beta theta C gamma2..gamma7
pHat = s4.pHat;

%% ============= FIGURE 1: bootstrap histograms, L0/beta/theta/C =============
f = figure('Visible','off','Position',[100 100 950 700]);
tl = tiledlayout(f, 2, 2, 'Padding','compact','TileSpacing','compact');
labels = {'L0 (ms)','\beta (ms/Mbps)','\theta (ms\cdotMbps)','C (Mbps)'};
for i = 1:4
    nexttile(tl);
    histogram(boot(:,i), 40, 'FaceColor',[0.16 0.47 0.84], 'EdgeColor','none');
    hold on;
    xline(pHat(i), '-', 'point estimate', 'Color',[0.11 0.69 0.48], 'LineWidth',2);
    xline(prctile(boot(:,i),2.5), ':', 'Color',[0.3 0.3 0.3]);
    xline(prctile(boot(:,i),97.5), ':', 'Color',[0.3 0.3 0.3]);
    title(labels{i}); grid on; box on;
end
title(tl, sprintf('Bootstrap distributions (B = %d) with the 95%% percentile interval', size(boot,1)));
saveas(f, fullfile(figDir, 'step04_bootstrap_histograms.png'));
close(f);

%% ============= FIGURE 2: gamma_q bootstrap intervals =============
gammaBootLo = [0, prctile(boot(:,5:10), 2.5, 1)];
gammaBootHi = [0, prctile(boot(:,5:10), 97.5, 1)];
gammaPoint  = [0, pHat(5:10)];

f = figure('Visible','off','Position',[100 100 900 550]);
x = 1:7;
errorbar(x, gammaPoint, gammaPoint-gammaBootLo, gammaBootHi-gammaPoint, 'o', ...
    'Color',[0.16 0.47 0.84], 'MarkerFaceColor',[0.16 0.47 0.84], 'LineWidth',1.8, 'CapSize',10);
set(gca,'XTick',1:7,'XTickLabel',arrayfun(@(k)sprintf('Q%d',k),1:7,'UniformOutput',false));
ylabel('\gamma_q (ms)');
title('Location offsets with bootstrap 95% intervals (Q1 = reference)');
grid on; box on;
saveas(f, fullfile(figDir, 'step04_gamma_bootstrap_ci.png'));
close(f);

%% ============= FIGURE 3: interval width, bootstrap vs asymptotic =============
ci_asym = s4.ci95_asymptotic;
bootWidth = prctile(boot,97.5,1) - prctile(boot,2.5,1);
asymWidth = (ci_asym(:,2) - ci_asym(:,1))';
ratio = asymWidth(3:4) ./ bootWidth(3:4);

f = figure('Visible','off','Position',[100 100 950 550]);
b1 = bar([bootWidth(3:4); asymWidth(3:4)]');
b1(1).FaceColor = [0.16 0.47 0.84];
b1(2).FaceColor = [0.92 0.41 0.20];
set(gca, 'YScale','log', 'XTickLabel', {'\theta (ms\cdotMbps)','C (Mbps)'});
ylabel('Width of the 95% interval (log scale)');
legend({'Bootstrap','Asymptotic (Jacobian)'}, 'Location','northwest');
title({'Width of the 95% intervals of \theta and C', ...
       sprintf('asymptotic / bootstrap width: \\theta %.1f, C %.1f', ratio(1), ratio(2))});
grid on; box on;
saveas(f, fullfile(figDir, 'step04_ci_width_comparison.png'));
close(f);

%% ============= FIGURE 4: theta vs C bootstrap scatter =============
r = corr(boot(:,3), boot(:,4));
f = figure('Visible','off','Position',[100 100 900 700]);
scatter(boot(:,4), boot(:,3), 14, [0.16 0.47 0.84], 'filled', 'MarkerFaceAlpha', 0.35);
hold on;
plot(pHat(4), pHat(3), 'p', 'MarkerSize', 16, 'MarkerFaceColor',[0.11 0.69 0.48], 'MarkerEdgeColor','k');
xlabel('C (Mbps)'); ylabel('\theta (ms\cdotMbps)');
title(sprintf('Bootstrap estimates of \\theta versus C (n = %d, r = %.2f)', size(boot,1), r));
legend({'Bootstrap replicates','Point estimate'}, 'Location','northwest');
grid on; box on;
saveas(f, fullfile(figDir, 'step04_theta_vs_C.png'));
close(f);

fprintf('Figures written to %s\n', figDir);
end
