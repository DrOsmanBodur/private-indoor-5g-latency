function step03_cross_validate()
% STEP 03 - Out-of-sample validation of the proposed model using the
% achieved throughput (Section 4.3.2, Table 6 row (a)).
%
% (a) Trial-based split: train = trials 1-3 (n=126), test = trials 4-5
%     (n=84). Fit on the training runs only and evaluate MAE/RMSE/R2 on the
%     test runs.
% (b) Leave-one-location-out (LOLO): for each of the 7 locations, fit on the
%     other 6 (with their own location offsets) and predict the held-out
%     location with the population-average offset of the 6 training
%     locations, since the offset of a new location is not known in advance.
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

opts = optimoptions('lsqcurvefit', 'Display','off', ...
    'MaxFunctionEvaluations', 2e5, 'MaxIterations', 2e4, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

pub.L0 = 8.64; pub.beta = 0.0305; pub.theta = 486.84; pub.C = 611.39;
pub.gamma = [0 -0.47 -1.45 -1.31 0.78 -0.48 -1.68];
pub.MAE = 2.52; pub.RMSE = 3.13; pub.R2 = 0.842;

%% ============= PART A: TRAIN/TEST SPLIT BY TRIAL =============
trainMask = trial <= 3;
testMask  = trial >= 4;
fprintf('\n=== Part A: Train (trial 1-3, n=%d) / Test (trial 4-5, n=%d) ===\n', ...
    sum(trainMask), sum(testMask));

[pA, gA] = fitQueueingModel(b(trainMask), qIdx(trainMask), L(trainMask), 7, opts);
trainPred = evalQueueingModel(pA, gA, b(trainMask), qIdx(trainMask));
testPred  = evalQueueingModel(pA, gA, b(testMask),  qIdx(testMask));

trainMetrics = computeMetrics(L(trainMask), trainPred);
testMetrics  = computeMetrics(L(testMask),  testPred);

fprintf('Fitted on TRAIN: L0=%.4f beta=%.4f theta=%.4f C=%.4f\n', pA(1),pA(2),pA(3),pA(4));
fprintf('  TRAIN (in-sample):     MAE=%.4f RMSE=%.4f R2=%.4f\n', trainMetrics.MAE, trainMetrics.RMSE, trainMetrics.R2);
fprintf('  TEST  (out-of-sample): MAE=%.4f RMSE=%.4f R2=%.4f\n', testMetrics.MAE, testMetrics.RMSE, testMetrics.R2);
fprintf('  (cf. in-sample fit on all 210 runs: MAE=%.4f RMSE=%.4f R2=%.4f)\n', 2.5205, 3.1330, 0.8422);
fprintf('  (cf. article, in-sample:            MAE=%.4f RMSE=%.4f R2=%.4f)\n', pub.MAE, pub.RMSE, pub.R2);

%% ============= PART B: LEAVE-ONE-LOCATION-OUT =============
Q = 7;
foldMAE = nan(Q,1); foldRMSE = nan(Q,1); foldR2 = nan(Q,1); foldN = nan(Q,1);
allResid = []; allPred = []; allActual = [];
foldL0 = nan(Q,1); foldBeta = nan(Q,1); foldTheta = nan(Q,1); foldC = nan(Q,1); foldMeanGammaOthers = nan(Q,1);

fprintf('\n=== Part B: Leave-One-Location-Out (LOLO) ===\n');
for e = 1:Q
    trainIdx = qIdx ~= e;
    testIdx  = qIdx == e;

    locsIncluded = setdiff(1:Q, e);
    qLocal = zeros(sum(trainIdx),1);
    qFullTrain = qIdx(trainIdx);
    for k = 1:numel(locsIncluded)
        qLocal(qFullTrain == locsIncluded(k)) = k;
    end

    [pE, gE] = fitQueueingModel(b(trainIdx), qLocal, L(trainIdx), numel(locsIncluded), opts);
    meanGammaOthers = mean(gE);

    bTest = b(testIdx); Ltest = L(testIdx);
    predTest = pE(1) + pE(2).*bTest + pE(3)./(pE(4) - bTest) + meanGammaOthers;

    m = computeMetrics(Ltest, predTest);
    foldMAE(e) = m.MAE; foldRMSE(e) = m.RMSE; foldR2(e) = m.R2; foldN(e) = sum(testIdx);
    foldL0(e) = pE(1); foldBeta(e) = pE(2); foldTheta(e) = pE(3); foldC(e) = pE(4);
    foldMeanGammaOthers(e) = meanGammaOthers;

    allResid  = [allResid;  Ltest - predTest]; %#ok<AGROW>
    allPred   = [allPred;   predTest]; %#ok<AGROW>
    allActual = [allActual; Ltest]; %#ok<AGROW>

    fprintf('  Q%d held out (n=%d): L0=%.3f beta=%.4f theta=%.2f C=%.2f | MAE=%.4f RMSE=%.4f R2=%.4f\n', ...
        e, foldN(e), pE(1), pE(2), pE(3), pE(4), m.MAE, m.RMSE, m.R2);
end

pooled = computeMetrics(allActual, allPred);
fprintf('\nPooled LOLO out-of-sample (all 7 folds, n=%d): MAE=%.4f RMSE=%.4f R2=%.4f\n', ...
    numel(allActual), pooled.MAE, pooled.RMSE, pooled.R2);

%% ============= SAVE OUTPUTS =============
trainTestTable = table({'Train(in-sample)';'Test(out-of-sample)'}, ...
    [trainMetrics.MAE; testMetrics.MAE], [trainMetrics.RMSE; testMetrics.RMSE], [trainMetrics.R2; testMetrics.R2], ...
    'VariableNames', {'Split','MAE','RMSE','R2'});
trainParamsTable = table({'L0';'beta';'theta';'C'}, pA(1:4)', ...
    'VariableNames', {'Param','Value'});
gammaTrainTable = table((1:7)', gA(:), 'VariableNames', {'qIdx','gamma_train_ms'});

loloPerFoldTable = table((1:Q)', foldN, foldL0, foldBeta, foldTheta, foldC, foldMeanGammaOthers, foldMAE, foldRMSE, foldR2, ...
    'VariableNames', {'heldOut_qIdx','n_test','L0','beta','theta','C','meanGammaOthers','MAE','RMSE','R2'});
loloPooledTable = table(pooled.MAE, pooled.RMSE, pooled.R2, numel(allActual), ...
    'VariableNames', {'MAE','RMSE','R2','n_total'});

writetable(trainTestTable, fullfile(here, 'stage3_trainTest_metrics.csv'));
writetable(trainParamsTable, fullfile(here, 'stage3_trainTest_params.csv'));
writetable(gammaTrainTable, fullfile(here, 'stage3_trainTest_gamma.csv'));
writetable(loloPerFoldTable, fullfile(here, 'stage3_lolo_perfold.csv'));
writetable(loloPooledTable, fullfile(here, 'stage3_lolo_pooled_metrics.csv'));

results = struct();
results.trainParams = pA; results.trainGamma = gA;
results.trainMetrics = trainMetrics; results.testMetrics = testMetrics;
results.loloPerFold = loloPerFoldTable; results.loloPooled = pooled;
results.n = n;
save(fullfile(here, 'stage3_results.mat'), 'results');

xlsFile = fullfile(here, 'stage3_results.xlsx');
if isfile(xlsFile), delete(xlsFile); end
writetable(trainTestTable, xlsFile, 'Sheet', 'TrainTest Metrics');
writetable(trainParamsTable, xlsFile, 'Sheet', 'TrainTest Params');
writetable(gammaTrainTable, xlsFile, 'Sheet', 'TrainTest Gamma');
writetable(loloPerFoldTable, xlsFile, 'Sheet', 'LOLO Per Fold');
writetable(loloPooledTable, xlsFile, 'Sheet', 'LOLO Pooled');

fprintf('\nSaved: stage3_trainTest_metrics.csv, stage3_trainTest_params.csv, stage3_trainTest_gamma.csv,\n');
fprintf('       stage3_lolo_perfold.csv, stage3_lolo_pooled_metrics.csv, stage3_results.mat, stage3_results.xlsx\n');
fprintf('Done.\n');

end

%% ============ LOCAL FUNCTIONS ============
function [pHat, gamma] = fitQueueingModel(b, qIdxLocal, L, nLoc, opts)
% Fits L(b,q) = L0 + beta*b + theta/(C-b) + gamma(q), gamma(1)=0 reference,
% for a generic set of nLoc location groups (qIdxLocal in 1..nLoc).
    modelFun = @(p, X) qModel(p, X(:,1), X(:,2), nLoc);
    nFree = nLoc - 1;
    theta0 = [8.64, 0.0305, 486.84, max(b)+1, zeros(1,nFree)];
    lb     = [0,    0,      0,      max(b)+1, -50*ones(1,nFree)];
    ub     = [50,   5,      5000,   5000,      50*ones(1,nFree)];
    X = [b, qIdxLocal];
    pHat = lsqcurvefit(modelFun, theta0, X, L, lb, ub, opts);
    gamma = [0, pHat(5:end)];
end

function y = qModel(p, b, qIdxLocal, nLoc)
    b = b(:); qIdxLocal = round(qIdxLocal(:));
    L0 = p(1); beta = p(2); theta = p(3); C = p(4);
    gamma = [0, p(5:3+nLoc)];
    qIdxLocal(qIdxLocal<1) = 1; qIdxLocal(qIdxLocal>nLoc) = nLoc;
    y = L0 + beta.*b + theta./(C - b) + gamma(qIdxLocal)';
end

function pred = evalQueueingModel(pHat, gamma, b, qIdxLocal)
    b = b(:); qIdxLocal = round(qIdxLocal(:));
    qIdxLocal(qIdxLocal<1) = 1; qIdxLocal(qIdxLocal>numel(gamma)) = numel(gamma);
    pred = pHat(1) + pHat(2).*b + pHat(3)./(pHat(4) - b) + gamma(qIdxLocal)';
end

function m = computeMetrics(actual, pred)
    res = actual - pred;
    m.MAE = mean(abs(res));
    m.RMSE = sqrt(mean(res.^2));
    m.R2 = 1 - sum(res.^2)/sum((actual - mean(actual)).^2);
end
