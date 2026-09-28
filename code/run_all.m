% RUN_ALL  Runs the MATLAB analysis pipeline (steps 01-10) and the MATLAB
% figure scripts in order, and prints a summary at the end.
%
% Inputs are read from ../data and ../raw. Outputs are written to this
% folder (steps) and to figures/ (figure scripts).
% Requires MATLAB R2025b with the Optimization Toolbox and the Statistics
% and Machine Learning Toolbox.

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir, fullfile(codeDir, 'figures'));

scripts = {
    'step01_parse_achieved_bandwidth'
    'step02_fit_model'
    'step03_cross_validate'
    'step04_bootstrap_ci'
    'step05_tail_jitter_stats'
    'step06_gamma_bandwidth_corr'
    'step07_metallic_flag'
    'step08_baseline_comparison'
    'step09_urllc_supplementary'
    'step10_imputation_sensitivity'
    'make_fig5_trial_latency'
    'make_fig6_nonlinear_term'
    'make_fig7_location_offsets'
    'make_fig8_measured_vs_predicted'
    'diag_step01_achieved_bandwidth'
    'diag_step02_fit'
    'diag_step03_cross_validation'
    'diag_step04_bootstrap'
    'diag_step05_jitter_loss'
    'diag_step06_gamma_bandwidth'
    };

runStatus = strings(numel(scripts), 1);
for k = 1:numel(scripts)
    fprintf('\n==================== %s ====================\n', scripts{k});
    t0 = tic;
    try
        feval(scripts{k});
        runStatus(k) = sprintf('OK     (%.1f s)', toc(t0));
    catch ME
        runStatus(k) = "FAILED: " + ME.message;
    end
end

fprintf('\n==================== SUMMARY ====================\n');
for k = 1:numel(scripts)
    fprintf('%-34s %s\n', scripts{k}, runStatus(k));
end
