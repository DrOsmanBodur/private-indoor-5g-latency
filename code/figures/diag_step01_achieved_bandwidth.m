function diag_step01_achieved_bandwidth()
% Diagnostic figures for step 01 (achieved throughput) and the TCP test.
%   1. achieved (received) versus target UDP throughput per target rate
%   2. deviation of the achieved from the target throughput (%)
%   3. TCP run-level mean throughput per location (Section 3.1: first
%      one-second interval excluded, mean of the remaining 59 intervals)
%      next to the UDP achieved throughput at the 500 Mbps target
% The TCP statistics of Section 3.1 are printed.
%
% Inputs:  data/achieved_bandwidth_210.csv, raw/raw_throughput_logs.zip
% Outputs: PNG files in figures/diagnostics/

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
figDir = fullfile(here, 'diagnostics');
if ~exist(figDir, 'dir'), mkdir(figDir); end

A = readtable(fullfile(root, 'data', 'achieved_bandwidth_210.csv'));
bwLevels = [1 10 50 100 200 500];
nB = numel(bwLevels);
Q = 7;

%% ============= FIGURE 1: achieved vs target throughput per target rate =============
f = figure('Visible','off','Position',[100 100 950 600]);
grp = A.B_target_Mbps;
boxplot(A.achieved_Mbps_server, grp, 'Labels', arrayfun(@(k)sprintf('%dM',k),bwLevels,'UniformOutput',false));
hold on;
plot(1:nB, bwLevels, 'rd', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName','Target');
set(gca, 'YScale', 'log');
ylabel('Throughput (Mbps, log scale)'); xlabel('Target UDP rate');
title('Achieved (received) versus target UDP throughput, n = 210');
legend('Location','northwest');
grid on; box on;
saveas(f, fullfile(figDir, 'step01_achieved_vs_target.png'));
close(f);

%% ============= FIGURE 2: deviation of achieved from target (%) =============
pctDev = 100*(A.achieved_Mbps_server - A.B_target_Mbps) ./ A.B_target_Mbps;
f = figure('Visible','off','Position',[100 100 900 550]);
boxplot(pctDev, grp, 'Labels', arrayfun(@(k)sprintf('%dM',k),bwLevels,'UniformOutput',false));
yline(0, 'k--');
ylabel('Achieved versus target deviation (%)'); xlabel('Target UDP rate');
title('Deviation of the achieved from the target UDP throughput');
grid on; box on;
saveas(f, fullfile(figDir, 'step01_achieved_target_deviation.png'));
close(f);

%% ============= FIGURE 3: TCP run-level means per location vs UDP at 500 Mbps =============
tcp = tcpRunMeans(fullfile(root, 'raw', 'raw_throughput_logs.zip'));   % 7 x 5
tcpMean = mean(tcp, 2)'; tcpMin = min(tcp, [], 2)'; tcpMax = max(tcp, [], 2)';

sel500 = A.B_target_Mbps == 500;
udpMean = arrayfun(@(qi) mean(A.achieved_Mbps_server(sel500 & A.qIdx==qi)), 1:Q);
udpMin  = arrayfun(@(qi) min(A.achieved_Mbps_server(sel500 & A.qIdx==qi)), 1:Q);
udpMax  = arrayfun(@(qi) max(A.achieved_Mbps_server(sel500 & A.qIdx==qi)), 1:Q);

f = figure('Visible','off','Position',[100 100 950 600]);
hold on;
errorbar((1:Q)-0.12, tcpMean, tcpMean-tcpMin, tcpMax-tcpMean, 'o', ...
    'Color',[0.75 0.1 0.1], 'MarkerFaceColor',[0.75 0.1 0.1], 'LineWidth',1.6, 'CapSize',8, ...
    'DisplayName','TCP test, run-level means (Figure 3)');
errorbar((1:Q)+0.12, udpMean, udpMean-udpMin, udpMax-udpMean, 's', ...
    'Color',[0.1 0.3 0.7], 'MarkerFaceColor',[0.1 0.3 0.7], 'LineWidth',1.6, 'CapSize',8, ...
    'DisplayName','UDP achieved throughput at the 500 Mbps target');
set(gca,'XTick',1:7,'XTickLabel',arrayfun(@(k)sprintf('Q%d',k),1:7,'UniformOutput',false));
ylabel('Throughput (Mbps)'); xlabel('Location');
title({'TCP throughput test and UDP runs at the 500 Mbps target (different tests)', ...
       'markers: mean over 5 runs; bars: minimum to maximum'});
legend('Location','southoutside'); grid on; box on;
saveas(f, fullfile(figDir, 'step01_tcp_vs_udp500_by_location.png'));
close(f);

%% ---- TCP statistics of Section 3.1 ----
fprintf('TCP run-level mean throughput (first interval excluded), n = %d runs\n', numel(tcp));
fprintf('  all runs: mean %.1f Mbps, range %.1f-%.1f Mbps\n', mean(tcp(:)), min(tcp(:)), max(tcp(:)));
for qi = 1:Q
    fprintf('  Q%d: mean %.1f  median %.1f  SD %.1f  IQR %.1f\n', qi, mean(tcp(qi,:)), ...
        median(tcp(qi,:)), std(tcp(qi,:)), iqr(tcp(qi,:)));
end
fprintf('Figures written to %s\n', figDir);
end

function M = tcpRunMeans(zipFile)
% Run-level mean TCP throughput per location (rows) and run (columns):
% mean of the one-second intervals after the first one.
tmp = tempname; unzip(zipFile, tmp);
cleanupTmp = onCleanup(@() rmdir(tmp, 's'));
M = nan(7, 5);
files = dir(fullfile(tmp, '**', 'k*_bandwidthResult*.txt'));
for i = 1:numel(files)
    tok = regexp(files(i).name, '^k(\d+)_bandwidthResult(\d+)\.txt$', 'tokens', 'once');
    fn = fullfile(files(i).folder, files(i).name);
    lines = splitlines(string(fileread(fn)));   % UTF-16 is detected from the byte-order mark
    t = regexp(lines, '^\[\s*\d+\]\s+([\d.]+)-([\d.]+)\s+sec\s+[\d.]+\s+\w?Bytes\s+([\d.]+)\s+(G|M|K)bits/sec', 'tokens', 'once');
    v = [];
    for k = 1:numel(t)
        if isempty(t{k}), continue; end
        s = str2double(t{k}(1)); e = str2double(t{k}(2));
        if e - s > 1.01, continue; end                     % skip the whole-run summary line
        scale = containers.Map({'G','M','K'}, {1000, 1, 1e-3});
        v(end+1, :) = [s, str2double(t{k}(3)) * scale(char(t{k}(4)))]; %#ok<AGROW>
    end
    v = sortrows(v, 1);
    M(str2double(tok{1}), str2double(tok{2})) = mean(v(2:end, 2));
end
end
