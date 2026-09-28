function step02_fit_model()
% STEP 02 - Full-data fit of the proposed model (Section 4.3.2, Table 5).
%   L(b,q) = L0 + beta*b + theta/(C-b) + gamma(q)
% b is the achieved UDP throughput reported by the iPerf server
% (achieved_Mbps_server), not the configured target rate.
%
% Also computes the residual norm with C held fixed at several values and
% the multi-start analysis over the initial value of C. The article values
% are printed for comparison. The output tables also keep the estimates of a
% preliminary fit with the target rate as predictor (columns *Stage0*); that
% fit is not used in the article.
%
% Reads data/latency_with_achieved_bw_210.csv; writes outputs to the code folder.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);

%% -------------------- LOAD DATA ------------------------
T = readtable(fullfile(root, 'data', 'latency_with_achieved_bw_210.csv'));

qIdx = T.qIdx(:);
b    = T.achieved_Mbps_server(:);   % achieved throughput, not the target rate
L    = T.latency_ms(:);

Q = 7;
n = numel(L);
fprintf('Loaded %d rows (expected 210).\n', n);
fprintf('Achieved bandwidth range: [%.2f, %.2f] Mbps\n', min(b), max(b));

%% -------------------- MODEL DEFINITION ------------------
modelFun = @(p, X) publishedModel(p, X(:,1), X(:,2));

theta0 = [8.64, 0.0305, 486.84, 611.39, zeros(1,6)];
lb     = [0,    0,      0,      max(b)+1, -50*ones(1,6)];
ub     = [50,   5,      5000,   5000,      50*ones(1,6)];

opts = optimoptions('lsqcurvefit', ...
    'Display','off', ...
    'MaxFunctionEvaluations', 2e5, ...
    'MaxIterations', 2e4, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12);

X = [b, qIdx];

if ~exist('lsqcurvefit', 'file')
    error('lsqcurvefit not found. Requires Optimization Toolbox.');
end

%% -------------------- FIT (IN-SAMPLE, all 210 points) ----
[pHat, resnorm, residual, exitflag, output] = lsqcurvefit( ...
    modelFun, theta0, X, L, lb, ub, opts);

fprintf('\nFit done. exitflag=%d, resnorm=%.4f\n', exitflag, resnorm);

L0    = pHat(1);
beta  = pHat(2);
theta = pHat(3);
C     = pHat(4);
gamma = [0, pHat(5:10)];

%% -------------------- METRICS -----------------------------
Lhat = modelFun(pHat, X);
res  = L - Lhat;

MAE  = mean(abs(res));
RMSE = sqrt(mean(res.^2));
SSres = sum(res.^2);
SStot = sum((L - mean(L)).^2);
R2   = 1 - SSres/SStot;

%% -------------------- IDENTIFIABILITY CHECK ----------------
% Refit with C held fixed at several values to show the theta-C ridge
% (flat residual norm).
Cgrid = [max(b)+1, 530, 550, 600, 611.39, 700, 800, 1000, 2000, 5000, 10000];
nC = numel(Cgrid);
profTheta = nan(nC,1); profResnorm = nan(nC,1);
for i = 1:nC
    Cfix = Cgrid(i);
    modelFunFixedC = @(p, X) publishedModelFixedC(p, X(:,1), X(:,2), Cfix);
    p0 = [8.64, 0.0305, 486.84, zeros(1,6)];
    lbF = [0, 0, 0, -50*ones(1,6)];
    ubF = [50, 5, 1e6, 50*ones(1,6)];
    [pF, rn] = lsqcurvefit(modelFunFixedC, p0, X, L, lbF, ubF, opts);
    profTheta(i) = pF(3);
    profResnorm(i) = rn;
end

%% -------------------- MULTI-START CHECK --------------------
C0list = [max(b)+1, 530, 611.39, 800, 1200, 2000, 5000];
nS = numel(C0list);
startC = nan(nS,1); hatC = nan(nS,1); hatTheta = nan(nS,1); hatL0 = nan(nS,1); hatBeta = nan(nS,1); hatResnorm = nan(nS,1); hatExit = nan(nS,1);
for i = 1:nS
    theta0i = [8.64, 0.0305, 486.84, C0list(i), zeros(1,6)];
    lbi = lb; lbi(4) = max(b)+1;
    [pI, rnI, ~, exI] = lsqcurvefit(modelFun, theta0i, X, L, lbi, ub, opts);
    startC(i) = C0list(i); hatC(i) = pI(4); hatTheta(i) = pI(3);
    hatL0(i) = pI(1); hatBeta(i) = pI(2); hatResnorm(i) = rnI; hatExit(i) = exI;
end

%% -------------------- REPORT -----------------------------
fprintf('\n=== This fit vs article vs preliminary fit with the target rate ===\n');
fprintf('%-10s %12s %12s %12s\n', 'Param', 'This fit', 'Article', 'Prelim.');
fprintf('%-10s %12.4f %12.4f %12.4f\n', 'L0',    L0,    8.64,   9.5067);
fprintf('%-10s %12.4f %12.4f %12.4f\n', 'beta',  beta,  0.0305, 0.0328);
fprintf('%-10s %12.4f %12.4f %12.4f\n', 'theta', theta, 486.84, 4.3948);
fprintf('%-10s %12.4f %12.4f %12.4f\n', 'C',     C,     611.39, 501.00);

fprintf('\n--- Location offsets gamma_q (Q1 fixed = 0 reference) ---\n');
gammaPub = [0,-0.47,-1.45,-1.31,0.78,-0.48,-1.68];
gammaS0  = [0,-0.4866,-1.4616,-1.3178,0.7174,-0.4884,-1.6980];
for qi = 1:Q
    fprintf('gamma_Q%d: this fit=%8.4f  article=%8.4f  prelim.=%8.4f\n', qi, gamma(qi), gammaPub(qi), gammaS0(qi));
end

fprintf('\n=== Metrics (in-sample, n=%d) ===\n', n);
fprintf('%-8s %12s %12s %12s\n', 'Metric', 'This fit', 'Article', 'Prelim.');
fprintf('%-8s %12.4f %12.4f %12.4f\n', 'MAE',  MAE,  2.52,   2.4996);
fprintf('%-8s %12.4f %12.4f %12.4f\n', 'RMSE', RMSE, 3.13,   3.1113);
fprintf('%-8s %12.4f %12.4f %12.4f\n', 'R2',   R2,   0.842,  0.8443);

fprintf('\n--- Multi-start check (does C still converge to one point?) ---\n');
fprintf('%-10s %10s %10s %10s %8s\n', 'C0_start', 'C_hat', 'theta_hat', 'resnorm', 'exit');
for i = 1:nS
    fprintf('%-10.2f %10.2f %10.4f %10.2f %8d\n', startC(i), hatC(i), hatTheta(i), hatResnorm(i), hatExit(i));
end

fprintf('\n--- C-fixed profile scan (identifiability ridge check) ---\n');
fprintf('%-10s %14s %10s\n', 'C_fixed', 'theta_at_C', 'resnorm');
for i = 1:nC
    fprintf('%-10.2f %14.4f %10.2f\n', Cgrid(i), profTheta(i), profResnorm(i));
end

%% -------------------- SAVE OUTPUTS -----------------------
paramTable = table({'L0';'beta';'theta';'C'}, [L0;beta;theta;C], [8.64;0.0305;486.84;611.39], [9.5067;0.0328;4.3948;501.00], ...
    'VariableNames', {'Param','Stage2_achieved_b','Published','Stage0_target_b'});
gammaTable = table((1:Q)', gamma(:), gammaPub(:), gammaS0(:), ...
    'VariableNames', {'qIdx','gamma_Stage2','gamma_Published','gamma_Stage0'});
metricsTable = table({'MAE';'RMSE';'R2'}, [MAE;RMSE;R2], [2.52;3.13;0.842], [2.4996;3.1113;0.8443], ...
    'VariableNames', {'Metric','Stage2_achieved_b','Published','Stage0_target_b'});
multiStartTable = table(startC, hatC, hatTheta, hatL0, hatBeta, hatResnorm, hatExit, ...
    'VariableNames', {'C0_start','C_hat','theta_hat','L0_hat','beta_hat','resnorm','exitflag'});
profileTable = table(Cgrid(:), profTheta, profResnorm, 'VariableNames', {'C_fixed','theta_at_C','resnorm'});

writetable(paramTable, fullfile(here, 'stage2_global_params.csv'));
writetable(gammaTable, fullfile(here, 'stage2_gamma_offsets.csv'));
writetable(multiStartTable, fullfile(here, 'stage2_multistart_check.csv'));
writetable(profileTable, fullfile(here, 'stage2_C_profile.csv'));

results = struct();
results.pHat = pHat; results.L0=L0; results.beta=beta; results.theta=theta; results.C=C;
results.gamma = gamma; results.MAE=MAE; results.RMSE=RMSE; results.R2=R2;
results.exitflag = exitflag; results.output = output; results.n = n;
results.bMin = min(b); results.bMax = max(b);
save(fullfile(here, 'stage2_results.mat'), 'results');

xlsFile = fullfile(here, 'stage2_results.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(metricsTable, xlsFile, 'Sheet', 'Metrics');
writetable(paramTable,   xlsFile, 'Sheet', 'Global Params');
writetable(gammaTable,   xlsFile, 'Sheet', 'Gamma Offsets');
writetable(multiStartTable, xlsFile, 'Sheet', 'MultiStart Check');
writetable(profileTable, xlsFile, 'Sheet', 'C Profile');

fprintf('\nSaved: stage2_global_params.csv, stage2_gamma_offsets.csv, stage2_multistart_check.csv,\n');
fprintf('       stage2_C_profile.csv, stage2_results.mat, stage2_results.xlsx\n');
fprintf('Done.\n');

end

%% ============ LOCAL FUNCTIONS ============
function y = publishedModel(p, b, qIdx)
    b = b(:); qIdx = round(qIdx(:));
    L0    = p(1);
    beta  = p(2);
    theta = p(3);
    C     = p(4);
    gammaFree = p(5:10);
    gamma = [0, gammaFree];
    qIdx(qIdx<1) = 1; qIdx(qIdx>7) = 7;
    y = L0 + beta.*b + theta./(C - b) + gamma(qIdx)';
end

function y = publishedModelFixedC(p, b, qIdx, Cfix)
    b = b(:); qIdx = round(qIdx(:));
    L0    = p(1);
    beta  = p(2);
    theta = p(3);
    gammaFree = p(4:9);
    gamma = [0, gammaFree];
    qIdx(qIdx<1) = 1; qIdx(qIdx>7) = 7;
    y = L0 + beta.*b + theta./(Cfix - b) + gamma(qIdx)';
end
