function make_fig8_measured_vs_predicted()
% Figure 8 of the article: measured and predicted mean latency per target rate.
%
% Measured: mean and standard deviation of the 35 runs per target rate.
% Predicted: mean of the full-data model predictions (same fit as
% step02_fit_model) over the same 35 runs. RMSE and R^2 are the in-sample
% metrics of the fit; MAXD is the largest difference between the rounded
% measured and predicted means. The values are checked against the
% published Figure 8.
%
% Input:   data/latency_with_achieved_bw_210.csv
% Output:  mean_measured_vs_predicted_latency.png (300 dpi) in this folder

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));

%% ---- full-data fit (Section 4.3.2) ----
T = readtable(fullfile(root, 'data', 'latency_with_achieved_bw_210.csv'));
b = T.achieved_Mbps_server(:); q = T.qIdx(:); L = T.latency_ms(:); B = T.B_Mbps(:);
modelFun = @(p, X) p(1) + p(2).*X(:,1) + p(3)./(p(4) - X(:,1)) + offsets(p, X(:,2));
opts = optimoptions('lsqcurvefit', 'Display','off', 'MaxFunctionEvaluations', 2e5, ...
    'MaxIterations', 2e4, 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
p0 = [8.64, 0.0305, 486.84, 611.39, zeros(1,6)];
lb = [0, 0, 0, max(b)+1, -50*ones(1,6)];
ub = [50, 5, 5000, 5000, 50*ones(1,6)];
p = lsqcurvefit(modelFun, p0, [b q], L, lb, ub, opts);
Lhat = modelFun(p, [b q]);

levels = [1 10 50 100 200 500]; n = numel(levels);
meas = zeros(1,n); sd = zeros(1,n); pred = zeros(1,n);
for k = 1:n
    m = B == levels(k);
    meas(k) = mean(L(m)); sd(k) = std(L(m)); pred(k) = mean(Lhat(m));
end
res = L - Lhat;
RMSE = sqrt(mean(res.^2));
R2 = 1 - sum(res.^2)/sum((L - mean(L)).^2);
MAXD = max(abs(round(meas, 2) - round(pred, 2)));

%% ---- figure ----
G = [21 194 134]/255; R = [236 15 71]/255;
fig = figure('Visible','off', 'Color','w', 'Units','pixels', 'Position', [50 50 1100 700]);
ax = axes(fig, 'Units','pixels', 'Position', [100 80 637.92 532.0]);
hold(ax, 'on');
x = 1:n;
h1 = errorbar(ax, x, meas, sd, '-o', 'Color', G, 'LineWidth', 2, 'MarkerSize', 7, ...
    'MarkerFaceColor', G, 'CapSize', 8);
h2 = plot(ax, x, pred, '--s', 'Color', R, 'LineWidth', 2.2, 'MarkerSize', 7.5, 'MarkerFaceColor', R);
xlim(ax, [0.7 6.3]); ylim(ax, [0 35]);
set(ax, 'XTick', x, 'XTickLabel', string(levels), 'FontSize', 11, 'LineWidth', 1.2, ...
    'TickDir','in', 'Box','off', 'XGrid','on', 'YGrid','on', 'GridLineStyle',':', ...
    'XColor','k', 'YColor','k', 'GridColor','k');
title(ax, 'Measured and Predicted Mean Latency Across Traffic Levels', 'FontWeight','bold');
xlabel(ax, 'Target UDP Bandwidth (Mbps)', 'FontWeight','bold');
ylabel(ax, 'Latency (ms)', 'FontWeight','bold');
legend(ax, [h1 h2], {'Measured mean \pm SD', 'Predicted mean'}, 'Location','northwest', ...
    'Box','off', 'FontSize', 9.9);

% side table
tx = {'FontName','Consolas', 'FontWeight','bold', 'FontSize',10, 'HorizontalAlignment','left', ...
      'VerticalAlignment','middle', 'Clipping','off', 'Interpreter','none'};
xT = 6.9367; yT = 32.387; dyT = 2.5317;
text(ax, xT, yT, ' BW | Measured | Predicted | SD', tx{:});
rows = strings(n, 1);
for k = 1:n
    rows(k) = sprintf('%3d | %8.2f | %9.2f | %5.2f', levels(k), meas(k), pred(k), sd(k));
    text(ax, xT, yT - k*dyT, rows(k), tx{:});
end
summary = [sprintf("RMSE | %.2f ms", RMSE); sprintf("R^2  | %.3f", R2); sprintf("MAXD | %.2f ms", MAXD)];
for k = 1:3
    text(ax, xT, 12.091 - (k-1)*1.848, summary(k), tx{:});
end

outFile = fullfile(here, 'mean_measured_vs_predicted_latency.png');
exportgraphics(fig, outFile, 'Resolution', 300);
close(fig);

%% ---- check against the published figure ----
published = [5.38 9.62 14.14 12.74 14.17 29.62; 8.81 9.11 10.46 12.15 15.60 29.54; 1.35 1.78 1.34 1.42 2.03 4.73];
same = isequal(round([meas; pred; sd], 2), published) && round(RMSE, 2) == 3.13 && ...
    round(R2, 3) == 0.842 && round(MAXD, 2) == 3.68;
fprintf('%s\n', rows);
fprintf('%s\n', summary);
fprintf('Same as the published Figure 8: %d\n', same);
fprintf('Saved %s\n', outFile);
end

function v = offsets(p, q)
    g = [0, p(5:10)];
    v = g(round(q))';
end
