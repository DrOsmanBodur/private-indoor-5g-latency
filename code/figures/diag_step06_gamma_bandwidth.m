function diag_step06_gamma_bandwidth()
% Diagnostic figure for step 06: location offsets gamma_q versus the
% location-wise variability of the achieved throughput (Section 5.2).
% Writes a self-contained HTML page with two inline SVG scatter plots.
%
% Inputs:  stage8_gamma_vs_bandwidthSD.csv and stage8_results.mat written by
%          step06_gamma_bandwidth_corr
% Output:  figures/diagnostics/step06_gamma_vs_bandwidth_sd.html

here = fileparts(mfilename('fullpath'));
codeDir = fileparts(here);
figDir = fullfile(here, 'diagnostics');
if ~exist(figDir, 'dir'), mkdir(figDir); end

T = readtable(fullfile(codeDir, 'stage8_gamma_vs_bandwidthSD.csv'));
load(fullfile(codeDir, 'stage8_results.mat'), 'results');

svgA = makeScatterSVG(T.bandwidth_SD_pooled_Mbps, T.gamma_q, T.qIdx, ...
    'Pooled throughput SD (Mbps)', 'gamma_q (ms)', ...
    sprintf('(a) Pooled over all target rates: r=%.3f, p=%.3f', results.r_pooled, results.p_pooled));

svgB = makeScatterSVG(T.bandwidth_SD_withinLevel_Mbps, T.gamma_q, T.qIdx, ...
    'Mean within-rate throughput SD (Mbps)', 'gamma_q (ms)', ...
    sprintf('(b) Within each target rate: r=%.3f, p=%.3f', results.r_within, results.p_within));

html = sprintf([...
'<!DOCTYPE html><html><head><meta charset="utf-8">' ...
'<title>gamma_q vs throughput SD</title>' ...
'<style>body{font-family:Arial,Helvetica,sans-serif;margin:24px;color:#222;}' ...
'.panels{display:flex;gap:24px;flex-wrap:wrap;} .panel{border:1px solid #ddd;border-radius:8px;padding:12px;}' ...
'h1{font-size:18px;} h2{font-size:14px;margin:4px 0 10px;} p{max-width:900px;}</style></head><body>' ...
'<h1>Location offset (gamma_q) versus achieved-throughput variability</h1>' ...
'<p>n = 7 locations. Panel (a) uses the SD over all 30 runs of a location, pooled over the 6 target rates; ' ...
'this measure is dominated by the spread between the target rates, which is the same at every location. ' ...
'Panel (b) uses the SD over the 5 trials at each target rate, averaged over the 6 rates, which isolates ' ...
'trial-to-trial throughput fluctuation at a fixed rate.</p>' ...
'<div class="panels">' ...
'<div class="panel">%s</div>' ...
'<div class="panel">%s</div>' ...
'</div></body></html>'], svgA, svgB);

fid = fopen(fullfile(figDir, 'step06_gamma_vs_bandwidth_sd.html'), 'w');
fwrite(fid, html, 'char');
fclose(fid);
fprintf('Saved %s\n', fullfile(figDir, 'step06_gamma_vs_bandwidth_sd.html'));
end

function svg = makeScatterSVG(x, y, labels, xlab, ylab, titleStr)
    W = 460; H = 380;
    mL = 60; mR = 20; mT = 30; mB = 55;
    plotW = W - mL - mR; plotH = H - mT - mB;

    xr = max(x) - min(x); if xr == 0, xr = 1; end
    yr = max(y) - min(y); if yr == 0, yr = 1; end
    xmin = min(x) - 0.12*xr; xmax = max(x) + 0.12*xr;
    ymin = min(y) - 0.12*yr; ymax = max(y) + 0.12*yr;

    toPxX = @(v) mL + (v - xmin) / (xmax - xmin) * plotW;
    toPxY = @(v) mT + (1 - (v - ymin) / (ymax - ymin)) * plotH;

    parts = {};
    parts{end+1} = sprintf('<svg viewBox="0 0 %d %d" width="100%%" style="max-width:%dpx" xmlns="http://www.w3.org/2000/svg">', W, H, W);
    parts{end+1} = sprintf('<text x="%d" y="16" font-size="12" font-weight="bold" fill="#222">%s</text>', mL, titleStr);
    % axes
    parts{end+1} = sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#888"/>', mL, mT+plotH, mL+plotW, mT+plotH);
    parts{end+1} = sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#888"/>', mL, mT, mL, mT+plotH);
    % zero line for gamma if within range
    if ymin < 0 && ymax > 0
        y0 = toPxY(0);
        parts{end+1} = sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#ccc" stroke-dasharray="3,3"/>', mL, y0, mL+plotW, y0);
    end
    % axis labels
    parts{end+1} = sprintf('<text x="%d" y="%d" font-size="11" fill="#444" text-anchor="middle">%s</text>', mL+plotW/2, H-8, xlab);
    parts{end+1} = sprintf('<text x="12" y="%.1f" font-size="11" fill="#444" text-anchor="middle" transform="rotate(-90 12 %.1f)">%s</text>', mT+plotH/2, mT+plotH/2, ylab);
    % tick labels (min/max)
    parts{end+1} = sprintf('<text x="%.1f" y="%d" font-size="9" fill="#666" text-anchor="start">%.2f</text>', mL, mT+plotH+16, min(x));
    parts{end+1} = sprintf('<text x="%.1f" y="%d" font-size="9" fill="#666" text-anchor="end">%.2f</text>', mL+plotW, mT+plotH+16, max(x));
    parts{end+1} = sprintf('<text x="%d" y="%.1f" font-size="9" fill="#666" text-anchor="end">%.2f</text>', mL-4, mT+4, max(y));
    parts{end+1} = sprintf('<text x="%d" y="%.1f" font-size="9" fill="#666" text-anchor="end">%.2f</text>', mL-4, mT+plotH, min(y));

    for i = 1:numel(x)
        px = toPxX(x(i)); py = toPxY(y(i));
        parts{end+1} = sprintf('<circle cx="%.1f" cy="%.1f" r="5" fill="#2563eb" fill-opacity="0.85"/>', px, py);
        parts{end+1} = sprintf('<text x="%.1f" y="%.1f" font-size="10" fill="#111">Q%d</text>', px+7, py-7, labels(i));
    end
    parts{end+1} = '</svg>';
    svg = strjoin(parts, '');
end
