function make_fig6_nonlinear_term()
% Figure 6 of the article: nonlinear contribution theta/(C-b) of the model.
%
% theta and C are estimated from the data with the full-data fit of
% step02_fit_model (Table 5); b_avg is the mean achieved throughput of the
% 35 runs per target rate. The values are compared with the published
% Figure 6.
%
% Input:   data/latency_with_achieved_bw_210.csv
% Output:  nonlinear_term_theta_over_C_minus_b.png (300 dpi) in this folder

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
theta = p(3); C = p(4);

levels = [1 10 50 100 200 500]; n = numel(levels);
bAvg = zeros(1,n);
for k = 1:n
    bAvg(k) = mean(b(B == levels(k)));
end
term = theta ./ (C - bAvg);

%% ---- figure ----
R = [232 23 74]/255; G = [29 191 133]/255;
W = 6.65; H = 5.54;                          % plot box, inches
ptx = 600 / (W*72); pty = 6 / (H*72);        % data units per point
fig = figure('Visible','off', 'Color','w', 'Units','inches', 'Position', [0.5 0.5 10.5 7]);
ax = axes(fig, 'Units','inches', 'Position', [0.8 0.7 W H]);
hold(ax, 'on');
bb = linspace(0, C - theta/6.3, 800);
plot(ax, bb, theta ./ (C - bb), '-', 'Color', R, 'LineWidth', 2.4);
plot(ax, bAvg, term, 'o', 'MarkerSize', 9.2, 'MarkerFaceColor', G, 'MarkerEdgeColor', G, ...
    'Clipping', 'off');
xlim(ax, [0 600]); ylim(ax, [0 6]);
set(ax, 'XTick', 0:100:600, 'YTick', 0:6, 'FontName','Arial', 'FontSize', 10.4, ...
    'LineWidth', 1, 'TickDir','in', 'TickLength', [4.8/(W*72) 4.8/(W*72)], 'Box','off', ...
    'XColor','k', 'YColor','k', 'XGrid','on', 'YGrid','on', 'GridLineStyle',':', ...
    'GridColor', [212 212 212]/255, 'GridAlpha', 1, 'Layer','bottom', ...
    'LabelFontSizeMultiplier', 11.5/10.4, 'TitleFontSizeMultiplier', 11.5/10.4);
xlabel(ax, 'Achieved Bandwidth, b (Mbps)', 'FontWeight','bold');
ylabel(ax, 'Nonlinear Contribution (ms)', 'FontWeight','bold');
title(ax, 'Nonlinear Contribution of the Term {\it\theta}/(C-b)', 'FontWeight','bold');

% point labels: offset from the marker (pt), horizontal and vertical alignment
labelPos = {[10 -22], 'left', 'top'; [6 24], 'left', 'bottom'; [18 -14], 'left', 'top';
            [14 20], 'center', 'bottom'; [8 -16], 'center', 'top'; [-12 0], 'right', 'middle'};
for k = 1:n
    d = labelPos{k,1};
    t = text(ax, bAvg(k) + d(1)*ptx, term(k) + d(2)*pty, sprintf('%d Mbps', levels(k)), ...
        'FontName','Arial', 'FontSize', 10, 'FontWeight','bold', 'BackgroundColor','w', ...
        'Margin', 1, 'HorizontalAlignment', labelPos{k,2}, 'VerticalAlignment', labelPos{k,3});
    e = t.Extent;                            % leader line from the label centre to the marker
    c = [(e(1) + e(3)/2)/ptx, (e(2) + e(4)/2)/pty];
    m = [bAvg(k)/ptx, term(k)/pty];
    u = (m - c) / norm(m - c);
    stop = m - 5*u;                          % ends 5 pt before the marker centre
    line(ax, [c(1) stop(1)]*ptx, [c(2) stop(2)]*pty, 'Color', [0.35 0.35 0.35], 'LineWidth', 0.8);
    uistack(t, 'top');
end

% side table: monospace columns starting 54 pt right of the plot box
tx = {'FontName','Consolas', 'FontWeight','bold', 'FontSize', 10.2, 'Interpreter','none', ...
      'HorizontalAlignment','left', 'VerticalAlignment','baseline', 'Clipping','off'};
adv = 0.5498 * 10.2;                         % Consolas character advance, pt
col = @(c) 600 + (54 + c*adv)*ptx;
row = @(y) 6 - y*pty;
yHead = 33; pitch = 28.9;
text(ax, col(0), row(yHead), ' BW |          |', tx{:});
text(ax, col(11), row(yHead), 'b', tx{:});
text(ax, col(12), row(yHead + 2.2), 'avg', tx{:}, 'FontSize', 0.7*10.2);
text(ax, col(19), row(yHead), [char(952) '/(C-b)'], tx{:});
rows = strings(n, 1);
for k = 1:n
    rows(k) = sprintf('%3d | %8.2f | %9.2f', levels(k), bAvg(k), term(k));
    text(ax, col(0), row(yHead + k*pitch), rows(k), tx{:});
end
params = [sprintf("%-4s | %.3f", char(952), theta); sprintf("%-4s | %.3f Mbps", 'C', C); ...
          sprintf("%-4s | %.2f ms", 'MAX', max(term))];
for k = 1:3
    text(ax, col(0), row(yHead + (n+2)*pitch + (k-1)*21.1), params(k), tx{:});
end

outFile = fullfile(here, 'nonlinear_term_theta_over_C_minus_b.png');
exportgraphics(fig, outFile, 'Resolution', 300);
close(fig);

%% ---- comparison with the published Figure 6 ----
% The published figure shows theta = 486.84 and C = 611.39 (as in the text
% of Section 5.1) and b_avg = 52.59 at 50 Mbps; the fit gives the Table 5
% values theta = 486.885 and C = 611.395, and the data give b_avg = 52.58.
pubBavg = [1.05 10.50 52.59 105.23 210.20 524.17];
pubTerm = [0.80 0.81 0.87 0.96 1.21 5.58];
fprintf('%s\n', rows);
fprintf('%s\n', strrep(params, char(952), 'theta'));
fprintf('theta/(C-b) and MAX same as the published Figure 6: %d\n', ...
    isequal(round(term, 2), pubTerm) && round(max(term), 2) == 5.58);
fprintf('theta and C same as Table 5: %d\n', round(theta, 3) == 486.885 && round(C, 3) == 611.395);
dif = find(round(bAvg, 2) ~= pubBavg);
for k = dif
    fprintf('b_avg at %d Mbps: %.2f (published Figure 6: %.2f)\n', levels(k), bAvg(k), pubBavg(k));
end
fprintf('Saved %s\n', outFile);
end

function v = offsets(p, q)
    g = [0, p(5:10)];
    v = g(round(q))';
end
