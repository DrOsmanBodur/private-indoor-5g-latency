function diag_step02_fit()
% Diagnostic figures for step 02 (full-data fit).
%   1. residual norm of the fit with C held fixed, over a range of C values
%      (the shallow theta-C ridge discussed in Section 5.1)
%   2. multi-start analysis: fitted C for each initial value C0 (Section 4.3.2)
%
% Inputs:  stage2_C_profile.csv, stage2_multistart_check.csv and
%          stage2_results.mat written by step02_fit_model
% Outputs: PNG files in figures/diagnostics/

here = fileparts(mfilename('fullpath'));
codeDir = fileparts(here);
figDir = fullfile(here, 'diagnostics');
if ~exist(figDir, 'dir'), mkdir(figDir); end

S2 = load(fullfile(codeDir, 'stage2_results.mat')); s2 = S2.results;
P2 = readtable(fullfile(codeDir, 'stage2_C_profile.csv'));
M2 = readtable(fullfile(codeDir, 'stage2_multistart_check.csv'));

%% ============= FIGURE 1: residual norm versus fixed C =============
f = figure('Visible','off','Position',[100 100 950 600]);
semilogx(P2.C_fixed, P2.resnorm, '-s', 'Color',[0.16 0.47 0.84], 'LineWidth',2, 'MarkerFaceColor',[0.16 0.47 0.84]);
xlabel('C held fixed (Mbps, log scale)'); ylabel('Residual norm (ms^2)');
title('Residual norm of the fit with C held fixed');
grid on; box on;
xline(s2.C, ':', sprintf('fitted C = %.2f Mbps', s2.C), 'LabelVerticalAlignment','bottom');
saveas(f, fullfile(figDir, 'step02_resnorm_vs_fixed_C.png'));
close(f);

%% ============= FIGURE 2: multi-start convergence =============
f = figure('Visible','off','Position',[100 100 950 550]);
semilogx(M2.C0_start, M2.C_hat, '-s', 'Color',[0.16 0.47 0.84], 'LineWidth',2, 'MarkerFaceColor',[0.16 0.47 0.84]);
xlabel('Initial value C_0 (Mbps, log scale)'); ylabel('Fitted C (Mbps)');
title('Multi-start analysis: fitted C for each initial value C_0');
grid on; box on;
ylim([480 630]);
saveas(f, fullfile(figDir, 'step02_multistart_convergence.png'));
close(f);

fprintf('Figures written to %s\n', figDir);
end
