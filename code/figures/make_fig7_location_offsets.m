function make_fig7_location_offsets()
% Figure 7 of the article: estimated location offsets gamma_q.
%
% The offsets come from the full-data fit of the proposed model (same
% settings as step02_fit_model). The values printed in the figure are
% checked against the values shown in the published Figure 7.
%
% Input:   data/latency_with_achieved_bw_210.csv
% Output:  location_offsets.png (300 dpi) in this folder

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));

%% ---- full-data fit (Section 4.3.2) ----
T = readtable(fullfile(root, 'data', 'latency_with_achieved_bw_210.csv'));
b = T.achieved_Mbps_server(:); q = T.qIdx(:); L = T.latency_ms(:);
modelFun = @(p, X) p(1) + p(2).*X(:,1) + p(3)./(p(4) - X(:,1)) + offsets(p, X(:,2));
opts = optimoptions('lsqcurvefit', 'Display','off', 'MaxFunctionEvaluations', 2e5, ...
    'MaxIterations', 2e4, 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
p0 = [8.64, 0.0305, 486.84, 611.39, zeros(1,6)];
lb = [0, 0, 0, max(b)+1, -50*ones(1,6)];
ub = [50, 5, 5000, 5000, 50*ones(1,6)];
p = lsqcurvefit(modelFun, p0, [b q], L, lb, ub, opts);
g = [0, p(5:10)];

%% ---- figure ----
C = [115 43 31; 204 51 41; 201 125 212; 110 92 242; 71 117 230; 36 191 224; 51 158 148] / 255;
fig = figure('Visible','off', 'Color','w', 'Units','pixels', 'Position', [50 50 1000 700]);
ax = axes(fig, 'Units','pixels', 'Position', [100 80 626.4 516.6]);
hold(ax, 'on');
for i = 1:7
    bar(ax, i, g(i), 0.55, 'FaceColor', C(i,:), 'FaceAlpha', 0.5, 'EdgeColor', C(i,:), 'LineWidth', 2);
end
yline(ax, 0, 'k', 'LineWidth', 1.6);
for i = 1:7
    s = 1; if g(i) < 0, s = -1; end
    text(ax, i, g(i) + s*0.13, sprintf('%+.2f', g(i)), 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', 'FontWeight','bold', 'FontSize', 10, 'Color', C(i,:));
end
xlim(ax, [0.4 7.6]);
ylim(ax, [min(g) - 0.45, max(g) + 0.45]);
set(ax, 'XTick', 1:7, 'XTickLabel', compose('Q%d', 1:7), 'FontSize', 11, 'LineWidth', 1.2, ...
    'TickDir','in', 'Box','off', 'XGrid','on', 'YGrid','on', 'GridLineStyle',':', ...
    'XColor','k', 'YColor','k', 'GridColor','k');
title(ax, 'Estimated Location-Dependent Latency Offsets', 'FontWeight','bold');
xlabel(ax, 'Measurement Points', 'FontWeight','bold');
ylabel(ax, 'Location Offset (ms)', 'FontWeight','bold');

% side table
tx = {'FontName','Consolas', 'FontWeight','bold', 'FontSize',10, 'HorizontalAlignment','left', ...
      'VerticalAlignment','middle', 'Clipping','off'};
xT = 8.4184; yT = 0.88957; dyT = 0.19916;
text(ax, xT, yT - 0.0073, ' Q  |  \gamma_q (ms)', tx{:});
rows = strings(7, 1);
for i = 1:7
    rows(i) = sprintf('Q%d  |  %+.2f', i, g(i));
    text(ax, xT, yT - i*dyT, rows(i), tx{:}, 'Color', C(i,:), 'Interpreter','none');
end
[gmax, imax] = max(g); [gmin, imin] = min(g);
summary = [sprintf("SPAN |  %.2f ms", gmax - gmin); sprintf("MAX  | Q%d  %+.2f", imax, gmax); ...
           sprintf("MIN  | Q%d  %+.2f", imin, gmin)];
for k = 1:3
    text(ax, xT, -0.97 - (k-1)*0.177, summary(k), tx{:}, 'Interpreter','none');
end

outFile = fullfile(here, 'location_offsets.png');
exportgraphics(fig, outFile, 'Resolution', 300);
close(fig);

%% ---- check against the published figure ----
published = ["+0.00" "-0.47" "-1.45" "-1.31" "+0.78" "-0.48" "-1.68"];
shown = compose('%+.2f', g);
fprintf('gamma_q shown: %s\n', strjoin(shown, '  '));
fprintf('Same as the published Figure 7: %d\n', isequal(string(shown), published) && ...
    summary(1) == "SPAN |  2.46 ms");
fprintf('Saved %s\n', outFile);
end

function v = offsets(p, q)
    g = [0, p(5:10)];
    v = g(round(q))';
end
