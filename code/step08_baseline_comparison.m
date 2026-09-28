function step08_baseline_comparison()
% STEP 08 - Out-of-sample comparison of the proposed model with an ablation
% model and three baselines (Section 5.4, Table 6).
%
%   (a)       proposed model:              L0 + beta*b + theta/(C-b) + gamma_q
%   (a-gamma) no-location ablation:        L0 + beta*b + theta/(C-b)
%   (b)       minimal M/M/1-type baseline: L0 + theta/(C-b)
%   (c)       polynomial regression poly22 in b and the location index
%   (d)       per-location saturation:     L0(q) + alpha(q) ln(1+b) + beta(q) max(0, b-Bth)
%
% (a) versus (a-gamma) isolates the contribution of the location offsets.
% All models are refitted from scratch with the achieved throughput, the
% same trial-based split and the same seven LOLO folds. The letters used
% here differ from Table 6 of the article: a = (a), a-gamma = (b),
% b = (c), c = (d), d = (e).
%
% Reads data/latency_with_achieved_bw_210.csv; writes outputs to the code folder.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);

T = readtable(fullfile(root, 'data', 'latency_with_achieved_bw_210.csv'));
qIdx  = T.qIdx(:);
trial = T.trial(:);
b     = T.achieved_Mbps_server(:);
L     = T.latency_ms(:);
n = numel(L);
fprintf('Loaded %d rows (expected 210).\n', n);
assert(n == 210, 'Expected 210 rows, got %d', n);
assert(~any(isnan(b)) && ~any(isnan(L)), 'NaN found in predictor or response.');

Bth = 200; % fixed threshold of the per-location saturation model (Mbps)
Q = 7;

opts = optimoptions('lsqcurvefit', 'Display','off', ...
    'MaxFunctionEvaluations', 2e5, 'MaxIterations', 2e4, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

models = {'a_proposed','ag_noLocation','b_minimalMM1','c_poly22','d_perLocationSat'};
nModels = numel(models);
nParamsMap = containers.Map(models, {10, 4, 3, 6, 21});

sanity = struct(); % accumulates all sanity-check results

%% ============= PART A: TRAIN/TEST SPLIT BY TRIAL =============
trainMask = trial <= 3;
testMask  = trial >= 4;
fprintf('\n=== Part A: Train (trial 1-3, n=%d) / Test (trial 4-5, n=%d) ===\n', ...
    sum(trainMask), sum(testMask));

trainTestResults = struct();
fittedParams = struct();
multistart = struct();

bTrain = b(trainMask); qTrain = qIdx(trainMask); LTrain = L(trainMask);
bTest  = b(testMask);  qTest  = qIdx(testMask);  LTest  = L(testMask);

% --- (a) Proposed full model ---
[pA, gA, resA] = fitPublished(bTrain, qTrain, LTrain, Q, opts);
predTr = evalPublished(pA, gA, bTrain, qTrain);
predTe = evalPublished(pA, gA, bTest,  qTest);
trainTestResults.a_proposed.train = computeMetrics(LTrain, predTr);
trainTestResults.a_proposed.test  = computeMetrics(LTest,  predTe);
fittedParams.a_proposed = struct('L0',pA(1),'beta',pA(2),'theta',pA(3),'C',pA(4),'gamma',gA);
checkNumeric(pA, 'a_proposed train params');
assert(pA(4) > max(bTrain), 'a_proposed: C must exceed max(b) in training data');
% multistart check
theta0_alt = [3, 0.08, 150, max(bTrain)+300, zeros(1,Q-1)];
lb = [0,0,0,max(bTrain)+1, -50*ones(1,Q-1)];
ub = [50,5,5000,5000, 50*ones(1,Q-1)];
modelFunA = @(p, X) qModel(p, X(:,1), X(:,2), Q);
[pA_alt, resA_alt] = lsqcurvefit(modelFunA, theta0_alt, [bTrain,qTrain], LTrain, lb, ub, opts);
multistart.a_proposed = struct('resnorm_main', resA, 'resnorm_alt', resA_alt, ...
    'maxParamDiff', max(abs(pA - pA_alt)), 'converged_same', max(abs(pA-pA_alt)) < 1e-2);

% --- (a-gamma) No-location ablation ---
[pAg, resAg] = fitNoLocation(bTrain, LTrain, opts);
predTr = evalNoLocation(pAg, bTrain);
predTe = evalNoLocation(pAg, bTest);
trainTestResults.ag_noLocation.train = computeMetrics(LTrain, predTr);
trainTestResults.ag_noLocation.test  = computeMetrics(LTest,  predTe);
fittedParams.ag_noLocation = struct('L0',pAg(1),'beta',pAg(2),'theta',pAg(3),'C',pAg(4));
checkNumeric(pAg, 'ag_noLocation train params');
assert(pAg(4) > max(bTrain), 'ag_noLocation: C must exceed max(b) in training data');
modelFunAg = @(p,x) p(1) + p(2).*x + p(3)./(p(4) - x);
theta0_alt2 = [3, 0.08, 150, max(bTrain)+300];
lb2 = [0,0,0,max(bTrain)+1]; ub2 = [50,5,5000,5000];
[pAg_alt, resAg_alt] = lsqcurvefit(modelFunAg, theta0_alt2, bTrain, LTrain, lb2, ub2, opts);
multistart.ag_noLocation = struct('resnorm_main', resAg, 'resnorm_alt', resAg_alt, ...
    'maxParamDiff', max(abs(pAg - pAg_alt)), 'converged_same', max(abs(pAg-pAg_alt)) < 1e-2);

% --- (b) Minimal M/M/1-type queueing baseline ---
[pB, resB] = fitMinimalMM1(bTrain, LTrain, opts);
predTr = evalMinimalMM1(pB, bTrain);
predTe = evalMinimalMM1(pB, bTest);
trainTestResults.b_minimalMM1.train = computeMetrics(LTrain, predTr);
trainTestResults.b_minimalMM1.test  = computeMetrics(LTest,  predTe);
fittedParams.b_minimalMM1 = struct('L0',pB(1),'theta',pB(2),'C',pB(3));
checkNumeric(pB, 'b_minimalMM1 train params');
assert(pB(3) > max(bTrain), 'b_minimalMM1: C must exceed max(b) in training data');
modelFunB = @(p,x) p(1) + p(2)./(p(3) - x);
theta0_alt3 = [3, 150, max(bTrain)+300];
lb3 = [0,0,max(bTrain)+1]; ub3 = [50,5000,5000];
[pB_alt, resB_alt] = lsqcurvefit(modelFunB, theta0_alt3, bTrain, LTrain, lb3, ub3, opts);
multistart.b_minimalMM1 = struct('resnorm_main', resB, 'resnorm_alt', resB_alt, ...
    'maxParamDiff', max(abs(pB - pB_alt)), 'converged_same', max(abs(pB-pB_alt)) < 1e-2);

% --- (c) Polynomial regression (poly22) - linear least squares, no multistart needed ---
pC = fitPoly22(bTrain, qTrain, LTrain);
predTr = evalPoly22(pC, bTrain, qTrain);
predTe = evalPoly22(pC, bTest,  qTest);
trainTestResults.c_poly22.train = computeMetrics(LTrain, predTr);
trainTestResults.c_poly22.test  = computeMetrics(LTest,  predTe);
fittedParams.c_poly22 = struct('coeffs', pC);
checkNumeric(pC, 'c_poly22 train params');

% --- (d) Per-location saturation model ---
[pD, resD] = fitPerLocation(bTrain, qTrain, LTrain, Q, Bth, opts);
predTr = evalPerLocation(pD, bTrain, qTrain, Q, Bth);
predTe = evalPerLocation(pD, bTest,  qTest,  Q, Bth);
trainTestResults.d_perLocationSat.train = computeMetrics(LTrain, predTr);
trainTestResults.d_perLocationSat.test  = computeMetrics(LTest,  predTe);
fittedParams.d_perLocationSat = struct('L0',pD(1:Q),'alpha',pD(Q+1:2*Q),'beta',pD(2*Q+1:3*Q));
checkNumeric(pD, 'd_perLocationSat train params');
modelFunD = @(p, X) perLocModel(p, X(:,1), X(:,2), Q, Bth);
theta0_alt4 = [2*ones(1,Q), 1*ones(1,Q), 0.1*ones(1,Q)];
lb4 = [zeros(1,Q), zeros(1,Q), zeros(1,Q)]; ub4 = [50*ones(1,Q), 20*ones(1,Q), 1*ones(1,Q)];
[pD_alt, resD_alt] = lsqcurvefit(modelFunD, theta0_alt4, [bTrain,qTrain], LTrain, lb4, ub4, opts);
multistart.d_perLocationSat = struct('resnorm_main', resD, 'resnorm_alt', resD_alt, ...
    'maxParamDiff', max(abs(pD - pD_alt)), 'converged_same', max(abs(pD-pD_alt)) < 1e-1);

for i = 1:nModels
    m = models{i};
    tr = trainTestResults.(m).train; te = trainTestResults.(m).test;
    fprintf('  %-16s (p=%2d): TRAIN MAE=%.4f RMSE=%.4f R2=%.4f | TEST MAE=%.4f RMSE=%.4f R2=%.4f\n', ...
        m, nParamsMap(m), tr.MAE, tr.RMSE, tr.R2, te.MAE, te.RMSE, te.R2);
end

%% ============= PART B: LEAVE-ONE-LOCATION-OUT =============
fprintf('\n=== Part B: Leave-One-Location-Out (LOLO) ===\n');

lolo = struct();
for i = 1:nModels
    lolo.(models{i}).allPred = [];
    lolo.(models{i}).allActual = [];
end
loloPerFold = struct(); % per-model, per-fold metrics
for i = 1:nModels
    loloPerFold.(models{i}) = repmat(struct('heldOutLoc',0,'n',0,'MAE',NaN,'RMSE',NaN,'R2',NaN), Q, 1);
end
spatialContribution = nan(Q,1); % mean gamma used for proposed model's held-out prediction, per fold

% Reference-location invariance sanity check accumulator (model a only)
refInvarianceMaxDiff = nan(Q,1);

for e = 1:Q
    trainIdx = qIdx ~= e;
    testIdx  = qIdx == e;
    bTr = b(trainIdx); qTr = qIdx(trainIdx); LTr = L(trainIdx);
    bTe = b(testIdx);  qTe = qIdx(testIdx);  LTe = L(testIdx);
    locsIncluded = setdiff(1:Q, e);
    nInc = numel(locsIncluded);

    % --- (a) Proposed: fit on 6 locations, predict held-out with population-average spatial effect ---
    qLocal = zeros(numel(qTr),1);
    for k = 1:nInc
        qLocal(qTr == locsIncluded(k)) = k;
    end
    [pE, gE] = fitPublished(bTr, qLocal, LTr, nInc, opts);
    meanGamma = mean(gE);
    spatialContribution(e) = meanGamma;
    predA = pE(1) + pE(2).*bTe + pE(3)./(pE(4) - bTe) + meanGamma;
    lolo.a_proposed.allPred = [lolo.a_proposed.allPred; predA];
    lolo.a_proposed.allActual = [lolo.a_proposed.allActual; LTe];
    fm = computeMetrics(LTe, predA);
    loloPerFold.a_proposed(e) = struct('heldOutLoc',e,'n',numel(LTe),'MAE',fm.MAE,'RMSE',fm.RMSE,'R2',fm.R2);

    % --- reference-location invariance sanity check for (a) ---
    refLoc1 = locsIncluded(1); refLoc2 = locsIncluded(end);
    predRef1 = fitPredictPopAvgWithRef(bTr, qTr, LTr, bTe, locsIncluded, refLoc1, opts);
    predRef2 = fitPredictPopAvgWithRef(bTr, qTr, LTr, bTe, locsIncluded, refLoc2, opts);
    refInvarianceMaxDiff(e) = max(abs(predRef1 - predRef2));

    % --- (a-gamma) no-location ablation: no location param -> generalizes directly ---
    pAg_e = fitNoLocation(bTr, LTr, opts);
    predAg = evalNoLocation(pAg_e, bTe);
    lolo.ag_noLocation.allPred = [lolo.ag_noLocation.allPred; predAg];
    lolo.ag_noLocation.allActual = [lolo.ag_noLocation.allActual; LTe];
    fm = computeMetrics(LTe, predAg);
    loloPerFold.ag_noLocation(e) = struct('heldOutLoc',e,'n',numel(LTe),'MAE',fm.MAE,'RMSE',fm.RMSE,'R2',fm.R2);

    % --- (b) minimal M/M/1-type: no location term -> generalizes directly ---
    pB_e = fitMinimalMM1(bTr, LTr, opts);
    predB = evalMinimalMM1(pB_e, bTe);
    lolo.b_minimalMM1.allPred = [lolo.b_minimalMM1.allPred; predB];
    lolo.b_minimalMM1.allActual = [lolo.b_minimalMM1.allActual; LTe];
    fm = computeMetrics(LTe, predB);
    loloPerFold.b_minimalMM1(e) = struct('heldOutLoc',e,'n',numel(LTe),'MAE',fm.MAE,'RMSE',fm.RMSE,'R2',fm.R2);

    % --- (c) poly22: q enters as continuous covariate -> use actual held-out qIdx ---
    pC_e = fitPoly22(bTr, qTr, LTr);
    predC = evalPoly22(pC_e, bTe, qTe);
    lolo.c_poly22.allPred = [lolo.c_poly22.allPred; predC];
    lolo.c_poly22.allActual = [lolo.c_poly22.allActual; LTe];
    fm = computeMetrics(LTe, predC);
    loloPerFold.c_poly22(e) = struct('heldOutLoc',e,'n',numel(LTe),'MAE',fm.MAE,'RMSE',fm.RMSE,'R2',fm.R2);

    % --- (d) per-location saturation: population-average (L0,alpha,beta) of remaining 6 ---
    pD_e = fitPerLocation(bTr, qLocal, LTr, nInc, Bth, opts);
    L0e = pD_e(1:nInc); alphaE = pD_e(nInc+1:2*nInc); betaE = pD_e(2*nInc+1:3*nInc);
    meanL0 = mean(L0e); meanAlpha = mean(alphaE); meanBeta = mean(betaE);
    predD = meanL0 + meanAlpha.*log(1+bTe) + meanBeta.*max(0, bTe-Bth);
    lolo.d_perLocationSat.allPred = [lolo.d_perLocationSat.allPred; predD];
    lolo.d_perLocationSat.allActual = [lolo.d_perLocationSat.allActual; LTe];
    fm = computeMetrics(LTe, predD);
    loloPerFold.d_perLocationSat(e) = struct('heldOutLoc',e,'n',numel(LTe),'MAE',fm.MAE,'RMSE',fm.RMSE,'R2',fm.R2);

    fprintf('  Q%d held out (n=%d) done.\n', e, sum(testIdx));
end

for i = 1:nModels
    m = models{i};
    pooled = computeMetrics(lolo.(m).allActual, lolo.(m).allPred);
    lolo.(m).pooled = pooled;
    fprintf('  %-16s LOLO pooled: MAE=%.4f RMSE=%.4f R2=%.4f\n', m, pooled.MAE, pooled.RMSE, pooled.R2);
end

sanity.refInvarianceMaxDiff_perFold = refInvarianceMaxDiff;
sanity.refInvarianceMaxDiff_overall = max(refInvarianceMaxDiff);
sanity.refInvariancePass = sanity.refInvarianceMaxDiff_overall < 0.05; % ms, nonlinear refit tolerance
fprintf('\nReference-location invariance check (model a, population-average prediction): max abs diff across folds = %.6g ms (pass<0.05: %d)\n', ...
    sanity.refInvarianceMaxDiff_overall, sanity.refInvariancePass);

%% ============= ASSEMBLE COMPARISON TABLE =============
labels = {'(a) Proposed full: L0+beta*b+theta/(C-b)+gamma_q'; ...
          '(a-gamma) No-location ablation: L0+beta*b+theta/(C-b)'; ...
          '(b) Minimal M/M/1-type queueing baseline: L0+theta/(C-b)'; ...
          '(c) Polynomial regression poly22(b,qIdx)'; ...
          '(d) Per-location saturation: L0(q)+alpha(q)ln(1+b)+beta(q)max(0,b-Bth)'};

nParams = zeros(nModels,1);
trMAE=zeros(nModels,1); trRMSE=zeros(nModels,1); trR2=zeros(nModels,1);
teMAE=zeros(nModels,1); teRMSE=zeros(nModels,1); teR2=zeros(nModels,1);
loMAE=zeros(nModels,1); loRMSE=zeros(nModels,1); loR2=zeros(nModels,1);

for i = 1:nModels
    m = models{i};
    nParams(i) = nParamsMap(m);
    trMAE(i)=trainTestResults.(m).train.MAE; trRMSE(i)=trainTestResults.(m).train.RMSE; trR2(i)=trainTestResults.(m).train.R2;
    teMAE(i)=trainTestResults.(m).test.MAE;  teRMSE(i)=trainTestResults.(m).test.RMSE;  teR2(i)=trainTestResults.(m).test.R2;
    loMAE(i)=lolo.(m).pooled.MAE; loRMSE(i)=lolo.(m).pooled.RMSE; loR2(i)=lolo.(m).pooled.R2;
end

compTable = table(labels, nParams, trMAE, trRMSE, trR2, teMAE, teRMSE, teR2, loMAE, loRMSE, loR2, ...
    'VariableNames', {'Model','nParams','train_MAE','train_RMSE','train_R2', ...
                       'test_MAE','test_RMSE','test_R2','LOLO_MAE','LOLO_RMSE','LOLO_R2'});
disp(' ');
disp(compTable);

%% ============= LOLO PER-FOLD TABLE =============
foldRows_model = {}; foldRows_loc = []; foldRows_n = []; foldRows_MAE = []; foldRows_RMSE = []; foldRows_R2 = [];
for i = 1:nModels
    m = models{i};
    for e = 1:Q
        r = loloPerFold.(m)(e);
        foldRows_model{end+1,1} = m; %#ok<AGROW>
        foldRows_loc(end+1,1) = r.heldOutLoc; %#ok<AGROW>
        foldRows_n(end+1,1) = r.n; %#ok<AGROW>
        foldRows_MAE(end+1,1) = r.MAE; %#ok<AGROW>
        foldRows_RMSE(end+1,1) = r.RMSE; %#ok<AGROW>
        foldRows_R2(end+1,1) = r.R2; %#ok<AGROW>
    end
end
loloPerFoldTable = table(foldRows_model, foldRows_loc, foldRows_n, foldRows_MAE, foldRows_RMSE, foldRows_R2, ...
    'VariableNames', {'Model','heldOutLoc','n','MAE','RMSE','R2'});

%% ============= LOCATION ABLATION EFFECT TABLE =============
teFull = trainTestResults.a_proposed.test; teAbl = trainTestResults.ag_noLocation.test;
loFull = lolo.a_proposed.pooled; loAbl = lolo.ag_noLocation.pooled;

ablation = table( ...
    {'test';'LOLO_pooled'}, ...
    [teAbl.MAE - teFull.MAE; loAbl.MAE - loFull.MAE], ...
    [teAbl.RMSE - teFull.RMSE; loAbl.RMSE - loFull.RMSE], ...
    [teFull.R2 - teAbl.R2; loFull.R2 - loAbl.R2], ...
    'VariableNames', {'Regime','DeltaMAE_ablation_minus_full','DeltaRMSE_ablation_minus_full','DeltaR2_full_minus_ablation'});
disp(' '); disp('Location-term (gamma_q) ablation effect (positive = full model better):'); disp(ablation);

% Check against the LOLO R2 of the proposed model reported in the article
expectedLOLO_R2 = 0.8175279;
actualLOLO_R2 = loR2(strcmp(models,'a_proposed'));
sanity.proposedLOLOR2_check = struct('expected', expectedLOLO_R2, 'actual', actualLOLO_R2, ...
    'absDiff', abs(actualLOLO_R2-expectedLOLO_R2), 'pass', abs(actualLOLO_R2-expectedLOLO_R2) < 1e-4);
fprintf('Proposed model LOLO R2 = %.7f (expected ~%.7f, article 0.818). Diff=%.2e. \n', ...
    actualLOLO_R2, expectedLOLO_R2, sanity.proposedLOLOR2_check.absDiff);

%% ============= SAVE OUTPUTS =============
writetable(compTable, fullfile(here, 'stage10_comparison_table.csv'));
writetable(loloPerFoldTable, fullfile(here, 'stage10_lolo_perfold.csv'));
writetable(ablation, fullfile(here, 'stage10_location_ablation_effect.csv'));

xlsFile = fullfile(here, 'stage10_results.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(compTable, xlsFile, 'Sheet', 'Comparison');
writetable(loloPerFoldTable, xlsFile, 'Sheet', 'LOLO_PerFold');
writetable(ablation, xlsFile, 'Sheet', 'Ablation_Effect');

paramRows_model = {}; paramRows_name = {}; paramRows_value = [];
pfn = fieldnames(fittedParams);
for i = 1:numel(pfn)
    m = pfn{i};
    s = fittedParams.(m);
    fn2 = fieldnames(s);
    for j = 1:numel(fn2)
        v = s.(fn2{j});
        for k = 1:numel(v)
            paramRows_model{end+1,1} = m; %#ok<AGROW>
            if numel(v) > 1
                paramRows_name{end+1,1} = sprintf('%s_%d', fn2{j}, k); %#ok<AGROW>
            else
                paramRows_name{end+1,1} = fn2{j}; %#ok<AGROW>
            end
            paramRows_value(end+1,1) = v(k); %#ok<AGROW>
        end
    end
end
paramTable = table(paramRows_model, paramRows_name, paramRows_value, ...
    'VariableNames', {'Model','Parameter','Value'});
writetable(paramTable, xlsFile, 'Sheet', 'Parameters');

sanityNames = {'proposed_LOLO_R2'; 'proposed_LOLO_R2_expected'; 'proposed_LOLO_R2_pass'; ...
    'refLocation_invariance_maxDiff_ms'; 'refLocation_invariance_pass'};
sanityValues = {actualLOLO_R2; expectedLOLO_R2; sanity.proposedLOLOR2_check.pass; ...
    sanity.refInvarianceMaxDiff_overall; sanity.refInvariancePass};
sanityTable = table(sanityNames, sanityValues, 'VariableNames', {'Check','Value'});
writetable(sanityTable, xlsFile, 'Sheet', 'Sanity_Checks', 'WriteVariableNames', true);

results = struct();
results.trainTestResults = trainTestResults;
results.lolo = lolo;
results.loloPerFold = loloPerFold;
results.compTable = compTable;
results.ablation = ablation;
results.fittedParams = fittedParams;
results.multistart = multistart;
results.sanity = sanity;
results.spatialContribution_perFold = spatialContribution;
save(fullfile(here, 'stage10_results.mat'), 'results');

fprintf('\nSaved: stage10_comparison_table.csv, stage10_lolo_perfold.csv, stage10_location_ablation_effect.csv, stage10_results.xlsx, stage10_results.mat\n');
fprintf('Done.\n');

end

%% ============ LOCAL FUNCTIONS ============

% ---- (a) Proposed full model: L0 + beta*b + theta/(C-b) + gamma(q) ----
function [pHat, gamma, resnorm] = fitPublished(b, qIdxLocal, L, nLoc, opts)
    modelFun = @(p, X) qModel(p, X(:,1), X(:,2), nLoc);
    nFree = nLoc - 1;
    theta0 = [8.64, 0.0305, 486.84, max(b)+1, zeros(1,nFree)];
    lb     = [0,    0,      0,      max(b)+1, -50*ones(1,nFree)];
    ub     = [50,   5,      5000,   5000,      50*ones(1,nFree)];
    X = [b, qIdxLocal];
    [pHat, resnorm] = lsqcurvefit(modelFun, theta0, X, L, lb, ub, opts);
    gamma = [0, pHat(5:end)];
end

function y = qModel(p, b, qIdxLocal, nLoc)
    b = b(:); qIdxLocal = round(qIdxLocal(:));
    L0 = p(1); beta = p(2); theta = p(3); C = p(4);
    gamma = [0, p(5:3+nLoc)];
    qIdxLocal(qIdxLocal<1) = 1; qIdxLocal(qIdxLocal>nLoc) = nLoc;
    y = L0 + beta.*b + theta./(C - b) + gamma(qIdxLocal)';
end

function pred = evalPublished(pHat, gamma, b, qIdxLocal)
    b = b(:); qIdxLocal = round(qIdxLocal(:));
    qIdxLocal(qIdxLocal<1) = 1; qIdxLocal(qIdxLocal>numel(gamma)) = numel(gamma);
    pred = pHat(1) + pHat(2).*b + pHat(3)./(pHat(4) - b) + gamma(qIdxLocal)';
end

% helper for reference-location invariance sanity check: fits the proposed
% model on a 6-location training set with a CHOSEN location mapped to the
% zero-offset reference, and returns the population-average prediction for
% the held-out test set. Used to verify the pop-average prediction does not
% depend on which of the 6 remaining locations is treated as reference.
function predAvg = fitPredictPopAvgWithRef(bTr, qOrigTr, LTr, bTe, locsIncluded, refLoc, opts)
    nInc = numel(locsIncluded);
    order = [refLoc, setdiff(locsIncluded, refLoc, 'stable')];
    qLocal = zeros(size(qOrigTr));
    for k = 1:nInc
        qLocal(qOrigTr == order(k)) = k;
    end
    [pE, gE] = fitPublished(bTr, qLocal, LTr, nInc, opts);
    meanGamma = mean(gE);
    predAvg = pE(1) + pE(2).*bTe + pE(3)./(pE(4) - bTe) + meanGamma;
end

% ---- (a-gamma) No-location ablation: L0 + beta*b + theta/(C-b), no gamma_q ----
function [pHat, resnorm] = fitNoLocation(b, L, opts)
    modelFun = @(p, x) p(1) + p(2).*x + p(3)./(p(4) - x);
    theta0 = [8.64, 0.0305, 486.84, max(b)+1];
    lb     = [0,    0,      0,      max(b)+1];
    ub     = [50,   5,      5000,   5000];
    [pHat, resnorm] = lsqcurvefit(modelFun, theta0, b, L, lb, ub, opts);
end

function pred = evalNoLocation(pHat, b)
    b = b(:);
    pred = pHat(1) + pHat(2).*b + pHat(3)./(pHat(4) - b);
end

% ---- (b) Minimal M/M/1-type queueing baseline: L0 + theta/(C-b), no linear or location term ----
function [pHat, resnorm] = fitMinimalMM1(b, L, opts)
    modelFun = @(p, x) p(1) + p(2)./(p(3) - x);
    theta0 = [8.64, 486.84, max(b)+1];
    lb     = [0,    0,      max(b)+1];
    ub     = [50,   5000,   5000];
    [pHat, resnorm] = lsqcurvefit(modelFun, theta0, b, L, lb, ub, opts);
end

function pred = evalMinimalMM1(pHat, b)
    b = b(:);
    pred = pHat(1) + pHat(2)./(pHat(3) - b);
end

% ---- (c) Polynomial regression baseline (poly22 in b and qIdx) - statistical, not physical ----
function pHat = fitPoly22(b, qIdxLocal, L)
    b = b(:); q = qIdxLocal(:); y = L(:);
    X = [ones(size(b)), b, q, b.^2, b.*q, q.^2];
    pHat = X \ y; % [p00; p10; p01; p20; p11; p02]
end

function pred = evalPoly22(pHat, b, qIdxLocal)
    b = b(:); q = qIdxLocal(:);
    X = [ones(size(b)), b, q, b.^2, b.*q, q.^2];
    pred = X * pHat;
end

% ---- (d) Per-location saturation model: L0(q) + alpha(q)*log(1+b) + beta(q)*max(0,b-Bth) ----
function [pHat, resnorm] = fitPerLocation(b, qIdxLocal, L, nLoc, Bth, opts)
    modelFun = @(p, X) perLocModel(p, X(:,1), X(:,2), nLoc, Bth);
    L0_0 = 4*ones(1,nLoc); alpha_0 = 1.8*ones(1,nLoc); beta_0 = 0.05*ones(1,nLoc);
    theta0 = [L0_0, alpha_0, beta_0];
    lb = [zeros(1,nLoc), zeros(1,nLoc), zeros(1,nLoc)];
    ub = [50*ones(1,nLoc), 20*ones(1,nLoc), 1*ones(1,nLoc)];
    X = [b, qIdxLocal];
    [pHat, resnorm] = lsqcurvefit(modelFun, theta0, X, L, lb, ub, opts);
end

function y = perLocModel(p, b, qIdxLocal, nLoc, Bth)
    b = b(:); q = round(qIdxLocal(:));
    q(q<1) = 1; q(q>nLoc) = nLoc;
    L0 = p(1:nLoc); alpha = p(nLoc+1:2*nLoc); beta = p(2*nLoc+1:3*nLoc);
    y = L0(q)' + alpha(q)'.*log(1+b) + beta(q)'.*max(0, b-Bth);
end

function pred = evalPerLocation(pHat, b, qIdxLocal, nLoc, Bth)
    pred = perLocModel(pHat, b, qIdxLocal, nLoc, Bth);
end

% ---- shared metric function ----
function m = computeMetrics(actual, pred)
    actual = actual(:); pred = pred(:);
    res = actual - pred;
    m.MAE = mean(abs(res));
    m.RMSE = sqrt(mean(res.^2));
    m.R2 = 1 - sum(res.^2)/sum((actual - mean(actual)).^2);
end

% ---- numeric sanity check on a fitted parameter vector ----
function checkNumeric(p, label)
    if any(isnan(p)) || any(isinf(p))
        error('Non-finite parameter detected in %s: %s', label, mat2str(p));
    end
end
