function step07_metallic_flag()
% STEP 07 - Metallic-environment indicator versus full location offsets
% (Section 5.2).
%
% Tests whether a single binary indicator for the locations in the area
% with metallic equipment and reflective surfaces (Q4, Q5, Q6) can replace
% the six free location offsets gamma_2..gamma_7, using the same trial-based
% split and leave-one-location-out (LOLO) scheme as steps 03 and 08, with
% the achieved throughput as predictor.
%
%   (a) Proposed model:     L = L0 + beta*b + theta/(C-b) + gamma_q            (10 params)
%   (e) Metallic indicator: L = L0 + beta*b + theta/(C-b) + delta*Metallic(q)  (5 params)
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

metallicSet = [4 5 6];
metallic = double(ismember(qIdx, metallicSet));
fprintf('Metallic(q): Q%s -> 1, others -> 0 (n_metallic=%d, n_other=%d)\n', ...
    mat2str(metallicSet), sum(metallic), sum(1-metallic));

opts = optimoptions('lsqcurvefit', 'Display','off', ...
    'MaxFunctionEvaluations', 2e5, 'MaxIterations', 2e4, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

models = {'a_published','e_metallicFlag'};
nModels = numel(models);

%% ============= PART A: TRAIN/TEST SPLIT BY TRIAL =============
trainMask = trial <= 3;
testMask  = trial >= 4;
fprintf('\n=== Part A: Train (trial 1-3, n=%d) / Test (trial 4-5, n=%d) ===\n', ...
    sum(trainMask), sum(testMask));

trainTestResults = struct();

% --- (a) Published queueing model ---
[pA, gA] = fitPublished(b(trainMask), qIdx(trainMask), L(trainMask), 7, opts);
predTr = evalPublished(pA, gA, b(trainMask), qIdx(trainMask));
predTe = evalPublished(pA, gA, b(testMask),  qIdx(testMask));
trainTestResults.a_published.nParams = 4 + (7-1);
trainTestResults.a_published.train = computeMetrics(L(trainMask), predTr);
trainTestResults.a_published.test  = computeMetrics(L(testMask),  predTe);

% --- (e) Metallic-flag variant ---
pE = fitMetallic(b(trainMask), metallic(trainMask), L(trainMask), opts);
predTr = evalMetallic(pE, b(trainMask), metallic(trainMask));
predTe = evalMetallic(pE, b(testMask),  metallic(testMask));
trainTestResults.e_metallicFlag.nParams = 5;
trainTestResults.e_metallicFlag.train = computeMetrics(L(trainMask), predTr);
trainTestResults.e_metallicFlag.test  = computeMetrics(L(testMask),  predTe);
trainTestResults.e_metallicFlag.delta = pE(5);

for i = 1:nModels
    m = models{i};
    tr = trainTestResults.(m).train; te = trainTestResults.(m).test;
    fprintf('  %-16s (p=%2d): TRAIN MAE=%.4f RMSE=%.4f R2=%.4f | TEST MAE=%.4f RMSE=%.4f R2=%.4f\n', ...
        m, trainTestResults.(m).nParams, tr.MAE, tr.RMSE, tr.R2, te.MAE, te.RMSE, te.R2);
end
fprintf('  delta (metallic coefficient, train fit) = %.4f ms\n', trainTestResults.e_metallicFlag.delta);

%% ============= PART B: LEAVE-ONE-LOCATION-OUT =============
Q = 7;
fprintf('\n=== Part B: Leave-One-Location-Out (LOLO) ===\n');

lolo = struct();
for i = 1:nModels
    lolo.(models{i}).allPred = [];
    lolo.(models{i}).allActual = [];
end
deltaPerFold = zeros(Q,1);

for eIdx = 1:Q
    trainIdx = qIdx ~= eIdx;
    testIdx  = qIdx == eIdx;
    bTr = b(trainIdx); qTr = qIdx(trainIdx); LTr = L(trainIdx);
    bTe = b(testIdx);  qTe = qIdx(testIdx);  LTe = L(testIdx);
    metTr = metallic(trainIdx); metTe = metallic(testIdx);

    locsIncluded = setdiff(1:Q, eIdx);
    qLocal = zeros(numel(qTr),1);
    for k = 1:numel(locsIncluded)
        qLocal(qTr == locsIncluded(k)) = k;
    end

    % (a) published: fit on 6 locations, predict held-out with mean gamma of the other 6
    [pAe, gAe] = fitPublished(bTr, qLocal, LTr, numel(locsIncluded), opts);
    meanGamma = mean(gAe);
    predA = pAe(1) + pAe(2).*bTe + pAe(3)./(pAe(4) - bTe) + meanGamma;
    lolo.a_published.allPred = [lolo.a_published.allPred; predA];
    lolo.a_published.allActual = [lolo.a_published.allActual; LTe];

    % (e) metallic flag: fit on 6 locations (their true Metallic values),
    % predict held-out location using ITS OWN true Metallic flag directly
    % (a genuinely new location's metallic-or-not status is knowable a
    % priori, unlike its individual gamma_q offset)
    pEe = fitMetallic(bTr, metTr, LTr, opts);
    predE = evalMetallic(pEe, bTe, metTe);
    lolo.e_metallicFlag.allPred = [lolo.e_metallicFlag.allPred; predE];
    lolo.e_metallicFlag.allActual = [lolo.e_metallicFlag.allActual; LTe];
    deltaPerFold(eIdx) = pEe(5);

    fprintf('  Q%d held out (n=%d) done.\n', eIdx, sum(testIdx));
end

for i = 1:nModels
    m = models{i};
    pooled = computeMetrics(lolo.(m).allActual, lolo.(m).allPred);
    lolo.(m).pooled = pooled;
    fprintf('  %-16s LOLO pooled: MAE=%.4f RMSE=%.4f R2=%.4f\n', m, pooled.MAE, pooled.RMSE, pooled.R2);
end
fprintf('  delta per LOLO fold: %s (mean=%.4f)\n', mat2str(round(deltaPerFold,3)), mean(deltaPerFold));

%% ============= ASSEMBLE SINGLE COMPARISON TABLE =============
labels = {'(a) Published queueing L0+beta*b+theta/(C-b)+gamma_q'; ...
          '(e) Metallic-flag L0+beta*b+theta/(C-b)+delta*Metallic(q)'};

nParams = zeros(nModels,1);
trMAE=zeros(nModels,1); trRMSE=zeros(nModels,1); trR2=zeros(nModels,1);
teMAE=zeros(nModels,1); teRMSE=zeros(nModels,1); teR2=zeros(nModels,1);
loMAE=zeros(nModels,1); loRMSE=zeros(nModels,1); loR2=zeros(nModels,1);

for i = 1:nModels
    m = models{i};
    nParams(i) = trainTestResults.(m).nParams;
    trMAE(i)=trainTestResults.(m).train.MAE; trRMSE(i)=trainTestResults.(m).train.RMSE; trR2(i)=trainTestResults.(m).train.R2;
    teMAE(i)=trainTestResults.(m).test.MAE;  teRMSE(i)=trainTestResults.(m).test.RMSE;  teR2(i)=trainTestResults.(m).test.R2;
    loMAE(i)=lolo.(m).pooled.MAE; loRMSE(i)=lolo.(m).pooled.RMSE; loR2(i)=lolo.(m).pooled.R2;
end

compTable = table(labels, nParams, trMAE, trRMSE, trR2, teMAE, teRMSE, teR2, loMAE, loRMSE, loR2, ...
    'VariableNames', {'Model','nParams','train_MAE','train_RMSE','train_R2', ...
                       'test_MAE','test_RMSE','test_R2','LOLO_MAE','LOLO_RMSE','LOLO_R2'});

disp(' ');
disp(compTable);

fprintf('\ntest R2 gap (a minus e): %.4f\n', teR2(1)-teR2(2));
fprintf('LOLO R2 gap (a minus e): %.4f\n', loR2(1)-loR2(2));

%% ============= SAVE OUTPUTS =============
writetable(compTable, fullfile(here, 'stage9_comparison_table.csv'));

xlsFile = fullfile(here, 'stage9_results.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(compTable, xlsFile, 'Sheet', 'Comparison');

results = struct();
results.trainTestResults = trainTestResults;
results.lolo = lolo;
results.compTable = compTable;
results.publishedParams = struct('L0',pA(1),'beta',pA(2),'theta',pA(3),'C',pA(4),'gamma',gA);
results.metallicParams = struct('L0',pE(1),'beta',pE(2),'theta',pE(3),'C',pE(4),'delta',pE(5));
results.deltaPerLOLOFold = deltaPerFold;
save(fullfile(here, 'stage9_results.mat'), 'results');

fprintf('\nSaved: stage9_comparison_table.csv, stage9_results.xlsx, stage9_results.mat\n');
fprintf('Done.\n');

end

%% ============ LOCAL FUNCTIONS ============

function [pHat, gamma] = fitPublished(b, qIdxLocal, L, nLoc, opts)
    p0 = [8, 0.03, 400, max(b)+50, zeros(1,nLoc-1)];
    lb = [-Inf, -Inf, 0, max(b)+1, -Inf(1,nLoc-1)];
    ub = [Inf, Inf, Inf, Inf, Inf(1,nLoc-1)];
    pHat = lsqcurvefit(@modelFun, p0, b, L, lb, ub, opts);
    gamma = [0, pHat(5:end)];

    function y = modelFun(p, bb)
        g = [0, p(5:3+nLoc)];
        y = p(1) + p(2).*bb + p(3)./(p(4) - bb) + g(qIdxLocal)';
    end
end

function pred = evalPublished(pHat, gamma, b, qIdxLocal)
    qIdxLocal(qIdxLocal<1) = 1; qIdxLocal(qIdxLocal>numel(gamma)) = numel(gamma);
    pred = pHat(1) + pHat(2).*b + pHat(3)./(pHat(4) - b) + gamma(qIdxLocal)';
end

function pHat = fitMetallic(b, metallic, L, opts)
    p0 = [8, 0.03, 400, max(b)+50, 0];
    lb = [-Inf, -Inf, 0, max(b)+1, -Inf];
    ub = [Inf, Inf, Inf, Inf, Inf];
    modelFun = @(p,bb) p(1) + p(2).*bb + p(3)./(p(4) - bb) + p(5).*metallic;
    pHat = lsqcurvefit(modelFun, p0, b, L, lb, ub, opts);
end

function pred = evalMetallic(pHat, b, metallic)
    pred = pHat(1) + pHat(2).*b + pHat(3)./(pHat(4) - b) + pHat(5).*metallic;
end

function m = computeMetrics(actual, pred)
    err = actual - pred;
    m.MAE  = mean(abs(err));
    m.RMSE = sqrt(mean(err.^2));
    ssRes = sum(err.^2);
    ssTot = sum((actual - mean(actual)).^2);
    m.R2 = 1 - ssRes/ssTot;
end
