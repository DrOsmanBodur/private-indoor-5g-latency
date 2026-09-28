function step01_parse_achieved_bandwidth()
% STEP 01 - Achieved throughput and Server Report fields from the raw iPerf
% UDP client logs.
%
% Parses the 210 raw iPerf logs in raw/raw_latency_logs.zip and extracts, per
% (location qIdx, target rate B_Mbps, trial):
%   - achieved_Mbps_client : client-side total-interval bandwidth (sent)
%   - achieved_Mbps_server : Server Report bandwidth (received); this is the
%                            achieved throughput used as the model predictor
%   - jitter_ms, loss_pct, latency avg/min/max/stdev from the Server Report line
%   - latency_corrupted_flag: one-way delay fields invalidated by the 32-bit
%                             counter wrap-around (Q1, 50 Mbps, trials 3-5)
%
% Does NOT modify any package file; all outputs go to the code folder.

outDir  = fileparts(mfilename('fullpath'));
rootDir = fileparts(outDir);
tmpDir  = tempname;
unzip(fullfile(rootDir, 'raw', 'raw_latency_logs.zip'), tmpDir);
cleanupTmp = onCleanup(@() rmdir(tmpDir, 's'));
dataDir = fullfile(tmpDir, 'latency');

files = dir(fullfile(dataDir, 'k*_latency*.txt'));
fprintf('Found %d raw files in %s\n', numel(files), dataDir);

nRows = numel(files);
qIdx        = zeros(nRows,1);
B_target    = zeros(nRows,1);
trial       = zeros(nRows,1);
ach_client  = nan(nRows,1);
ach_server  = nan(nRows,1);
jitter_ms   = nan(nRows,1);
loss_pct    = nan(nRows,1);
lat_avg     = nan(nRows,1);
lat_min     = nan(nRows,1);
lat_max     = nan(nRows,1);
lat_stdev   = nan(nRows,1);
corrupted   = false(nRows,1);
fname_col   = strings(nRows,1);

for i = 1:nRows
    fn = files(i).name;
    fname_col(i) = string(fn);

    % --- filename -> qIdx, B_target, trial ---
    tok = regexp(fn, '^k(\d+)_latency_(\d+)M_(\d+)\.txt$', 'tokens', 'once');
    if ~isempty(tok)
        qIdx(i)     = str2double(tok{1});
        B_target(i) = str2double(tok{2});
        trial(i)    = str2double(tok{3});
    else
        tok = regexp(fn, '^k(\d+)_latencyResult(\d+)\.txt$', 'tokens', 'once');
        if isempty(tok)
            warning('Unrecognized filename pattern, skipping: %s', fn);
            continue
        end
        qIdx(i)     = str2double(tok{1});
        B_target(i) = 1;   % "Result" files = 1 Mbps target (confirmed vs latency_210.csv)
        trial(i)    = str2double(tok{2});
    end

    % --- read file, auto-detecting encoding from BOM ---
    % Most files are UTF-16LE (with FF FE BOM), but a subset (all k4_latency_*M_*.txt,
    % i.e. location 4's bandwidth-limited trials) are plain ASCII/UTF-8 with no BOM.
    fid = fopen(fullfile(dataDir, fn), 'rb');
    bom2 = fread(fid, 2, 'uint8=>uint8')';
    fclose(fid);
    if isequal(bom2, uint8([255 254]))
        enc = 'UTF-16LE';
    elseif isequal(bom2, uint8([254 255]))
        enc = 'UTF-16BE';
    else
        enc = 'UTF-8';
    end
    if startsWith(enc, 'UTF-16')
        raw = fileread(fullfile(dataDir, fn));   % UTF-16 is detected from the byte-order mark
    else
        raw = fileread(fullfile(dataDir, fn), 'Encoding', enc);
    end
    lines = regexp(raw, '\r\n|\n|\r', 'split');

    % --- client total line: the line immediately before "[  1] Sent ... datagrams" ---
    sentIdx = find(~cellfun(@isempty, regexp(lines, '^\[\s*\d+\]\s+Sent\s+\d+\s+datagrams', 'once')), 1);
    if ~isempty(sentIdx) && sentIdx > 1
        totalLine = lines{sentIdx - 1};
        m = regexp(totalLine, '([\d.]+)\s*(Mbits|Kbits|Gbits)/sec', 'tokens', 'once');
        if ~isempty(m)
            val = str2double(m{1});
            unit = m{2};
            ach_client(i) = convertToMbps(val, unit);
        end
    end

    % --- server report data line: 2 lines after "Server Report:" ---
    srIdx = find(~cellfun(@isempty, regexp(lines, 'Server Report:', 'once')), 1);
    if ~isempty(srIdx) && (srIdx + 2) <= numel(lines)
        dataLine = lines{srIdx + 2};

        m = regexp(dataLine, '([\d.]+)\s*(Mbits|Kbits|Gbits)/sec', 'tokens', 'once');
        if ~isempty(m)
            ach_server(i) = convertToMbps(str2double(m{1}), m{2});
        end

        m = regexp(dataLine, '([\d.]+)\s*ms\s+(\d+)/(\d+)\s*\(([\d.]+)%\)', 'tokens', 'once');
        if ~isempty(m)
            jitter_ms(i) = str2double(m{1});
            loss_pct(i)  = str2double(m{4});
        end

        m = regexp(dataLine, '([\-\d.]+)/([\-\d.]+)/([\-\d.]+)/([\-\d.]+)\s*ms\s+\d+\s*pps', 'tokens', 'once');
        if ~isempty(m)
            lat_avg(i)   = str2double(m{1});
            lat_min(i)   = str2double(m{2});
            lat_max(i)   = str2double(m{3});
            lat_stdev(i) = str2double(m{4});
            % overflow/underflow bug: iperf reports a huge (~4.29e6 ms) value
            % when avg/min latency underflows (negative wraps to uint32)
            if abs(lat_avg(i)) > 1e5 || abs(lat_min(i)) > 1e5 || abs(lat_max(i)) > 1e5
                corrupted(i) = true;
            end
        end
    end
end

T = table(qIdx, B_target, trial, ach_client, ach_server, jitter_ms, loss_pct, ...
          lat_avg, lat_min, lat_max, lat_stdev, corrupted, fname_col, ...
          'VariableNames', {'qIdx','B_target_Mbps','trial','achieved_Mbps_client', ...
          'achieved_Mbps_server','jitter_ms_raw','loss_pct_raw','latency_avg_ms_raw', ...
          'latency_min_ms_raw','latency_max_ms_raw','latency_stdev_ms_raw', ...
          'latency_corrupted_flag','source_file'});
T = sortrows(T, {'qIdx','B_target_Mbps','trial'});

outCsv = fullfile(outDir, 'achieved_bandwidth_210.csv');
writetable(T, outCsv);
fprintf('Wrote %s (%d rows)\n', outCsv, height(T));

% --- sanity report ---
fprintf('\n--- Missing/failed parses ---\n');
fprintf('achieved_Mbps_client NaN count: %d\n', sum(isnan(T.achieved_Mbps_client)));
fprintf('achieved_Mbps_server NaN count: %d\n', sum(isnan(T.achieved_Mbps_server)));
fprintf('\n--- Corrupted latency rows (overflow bug in Server Report) ---\n');
disp(T(T.latency_corrupted_flag, {'qIdx','B_target_Mbps','trial','latency_avg_ms_raw','latency_min_ms_raw','source_file'}));

fprintf('\n--- Rows with NaN achieved bandwidth (parse failures) ---\n');
disp(T(isnan(T.achieved_Mbps_server), {'qIdx','B_target_Mbps','trial','source_file'}));

% --- merge with existing target-bandwidth latency_210.csv for downstream stages ---
lat210Path = fullfile(rootDir, 'data', 'latency_210.csv');
if isfile(lat210Path)
    L = readtable(lat210Path);
    Tjoin = renamevars(T, 'B_target_Mbps', 'B_Mbps');
    M = outerjoin(L, Tjoin, 'Keys', {'qIdx','B_Mbps','trial'}, 'RightVariables', ...
        {'achieved_Mbps_client','achieved_Mbps_server', ...
         'jitter_ms_raw','loss_pct_raw','latency_corrupted_flag'}, ...
        'MergeKeys', true);
    fprintf('\nRows after join (expect 210): %d\n', height(M));
    M = sortrows(M, {'qIdx','B_Mbps','trial'});
    mergedCsv = fullfile(outDir, 'latency_with_achieved_bw_210.csv');
    writetable(M, mergedCsv);
    fprintf('Wrote %s (%d rows)\n', mergedCsv, height(M));
end

end

function mbps = convertToMbps(val, unit)
switch unit
    case 'Gbits'
        mbps = val * 1000;
    case 'Mbits'
        mbps = val;
    case 'Kbits'
        mbps = val / 1000;
    otherwise
        mbps = NaN;
end
end
