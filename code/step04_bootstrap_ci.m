function step04_bootstrap_ci()
% STEP 04 - Parameter confidence intervals of the proposed model (Table 5),
% using the achieved throughput.
%
%   L(b,q) = L0 + beta*b + theta/(C-b) + gamma_q
%
% Two approaches, computed side by side:
%   (a) Nonparametric cell bootstrap: within each of the 42 (location,
%       target rate) cells the 5 trials are resampled with replacement
%       (preserving the balanced design), the model is refitted to each of
%       B = 1000 resampled datasets, and percentile intervals are taken.
%   (b) Asymptotic (Jacobian-based) intervals via nlparci, from the fit to
%       all 210 runs. These assume local linearity and are less reliable for
%       theta and C because of the shallow theta-C ridge (see step 02).
%
% Reads data/latency_with_achieved_bw_210.csv; writes outputs to the code folder.

clear; clc; rng(42); % reproducible bootstrap
here = fileparts(mfilename('fullpath'));
root = fileparts(here);

T = readtable(fullfile(root, 'data', 'latency_with_achieved_bw_210.csv'));
qIdx = T.qIdx(:); B_target = T.B_Mbps(:); trial = T.trial(:);
b = T.achieved_Mbps_server(:); L = T.latency_ms(:);
n = numel(L);
fprintf('Loaded %d rows (expected 210).\n', n);

opts = optimoptions('lsqcurvefit', 'Display','off', ...
    'MaxFunctionEvaluations', 2e5, 'MaxIterations', 2e4, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

pub.L0 = 8.64; pub.beta = 0.0305; pub.theta = 486.84; pub.C = 611.39;
pub.gamma = [0 -0.47 -1.45 -1.31 0.78 -0.48 -1.68];

%% ============= FULL-DATA FIT (point estimate + Jacobian for asymptotic CI) =============
modelFun = @(p, X) qModel(p, X(:,1), X(:,2));
theta0 = [8.64, 0.0305, 486.84, 611.39, zeros(1,6)];
lb = [0, 0, 0, max(b)+1, -50*ones(1,6)];
ub = [50, 5, 5000, 5000, 50*ones(1,6)];
X = [b, qIdx];

[pHat, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(modelFun, theta0, X, L, lb, ub, opts); %#ok<ASGLU>
fprintf('Full-data fit: L0=%.4f beta=%.4f theta=%.4f C=%.4f (resnorm=%.2f)\n', pHat(1),pHat(2),pHat(3),pHat(4), resnorm);

ci95 = nlparci(pHat, residual, 'jacobian', jacobian); % 10x2, columns [lo hi]

%% ============= CELL BOOTSTRAP =============
cells = unique([qIdx, B_target], 'rows');
nCells = size(cells,1);
cellIdx = cell(nCells,1);
for c = 1:nCells
    cellIdx{c} = find(qIdx==cells(c,1) & B_target==cells(c,2));
end

Bn = 1000;
bootParams = nan(Bn, 10); % L0 beta theta C gamma2..gamma7
bootExit = nan(Bn,1);
fprintf('Running %d cell-bootstrap replicates...\n', Bn);
tic;
for rep = 1:Bn
    idx = zeros(n,1);
    pos = 1;
    for c = 1:nCells
        cIdx = cellIdx{c};
        nC = numel(cIdx);
        samp = cIdx(randi(nC, nC, 1));
        idx(pos:pos+nC-1) = samp;
        pos = pos + nC;
    end
    bb = b(idx); qq = qIdx(idx); LL = L(idx);
    Xb = [bb, qq];
    try
        [pB, ~, ~, exB] = lsqcurvefit(modelFun, pHat, Xb, LL, lb, ub, opts);
        bootParams(rep,:) = pB;
        bootExit(rep) = exB;
    catch
        bootExit(rep) = -99;
    end
    if mod(rep,200)==0
        fprintf('  ...%d/%d done (%.1fs elapsed)\n', rep, Bn, toc);
    end
end
fprintf('Bootstrap done in %.1fs.\n', toc);

validReps = bootExit > 0;
fprintf('Valid replicates: %d / %d\n', sum(validReps), Bn);
bootParams = bootParams(validReps, :);

paramNames = {'L0','beta','theta','C','gamma2','gamma3','gamma4','gamma5','gamma6','gamma7'};
bootMean = mean(bootParams, 1);
bootSE   = std(bootParams, 0, 1);
bootCIlo = prctile(bootParams, 2.5, 1);
bootCIhi = prctile(bootParams, 97.5, 1);

%% ============= REPORT =============
fprintf('\n=== Parameter estimates with 95%% CIs (point estimate = full-data fit) ===\n');
fprintf('%-8s %10s | %20s | %20s\n', 'Param', 'PointEst', 'Bootstrap 95%% CI', 'Asymptotic (Jacobian) CI');
for i = 1:10
    fprintf('%-8s %10.4f | [%8.4f, %8.4f] | [%8.4f, %8.4f]\n', paramNames{i}, pHat(i), bootCIlo(i), bootCIhi(i), ci95(i,1), ci95(i,2));
end

fprintf('\n--- CI width comparison (bootstrap vs asymptotic) ---\n');
for i = 1:10
    wB = bootCIhi(i) - bootCIlo(i);
    wA = ci95(i,2) - ci95(i,1);
    fprintf('%-8s bootstrap width=%10.4f  asymptotic width=%10.4f  ratio(boot/asym)=%.2f\n', paramNames{i}, wB, wA, wB/wA);
end

%% ============= SAVE OUTPUTS =============
pointTable = table(paramNames', pHat(:), bootMean(:), bootSE(:), bootCIlo(:), bootCIhi(:), ci95(:,1), ci95(:,2), ...
    'VariableNames', {'Param','PointEstimate','BootstrapMean','BootstrapSE','Boot_CI_lo','Boot_CI_hi','Asymp_CI_lo','Asymp_CI_hi'});
writetable(pointTable, fullfile(here, 'stage4_parameter_CIs.csv'));

bootRawTable = array2table(bootParams, 'VariableNames', paramNames);
writetable(bootRawTable, fullfile(here, 'stage4_bootstrap_raw.csv'));

results = struct();
results.pHat = pHat; results.resnorm = resnorm; results.ci95_asymptotic = ci95;
results.bootParams = bootParams; results.bootMean = bootMean; results.bootSE = bootSE;
results.bootCIlo = bootCIlo; results.bootCIhi = bootCIhi; results.paramNames = {paramNames};
results.Bn = Bn; results.validReps = sum(validReps); results.n = n;
save(fullfile(here, 'stage4_results.mat'), 'results');

xlsFile = fullfile(here, 'stage4_results.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(pointTable, xlsFile, 'Sheet', 'Parameter CIs');

fprintf('\nSaved: stage4_parameter_CIs.csv, stage4_bootstrap_raw.csv, stage4_results.mat, stage4_results.xlsx\n');
fprintf('Done.\n');

end

%% ============ LOCAL FUNCTIONS ============
function y = qModel(p, b, qIdx)
    b = b(:); qIdx = round(qIdx(:));
    L0 = p(1); beta = p(2); theta = p(3); C = p(4);
    gamma = [0, p(5:10)];
    qIdx(qIdx<1) = 1; qIdx(qIdx>7) = 7;
    y = L0 + beta.*b + theta./(C - b) + gamma(qIdx)';
end
