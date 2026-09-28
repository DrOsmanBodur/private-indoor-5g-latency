function model = fit_latency_model(location, throughput, latency, opts)
% FIT_LATENCY_MODEL  Fits the latency model of the article to measurements.
%
%   L(b,q) = L0 + beta*b + theta/(C-b) + gamma_q,  gamma = 0 at the reference location
%
%   model = fit_latency_model(location, throughput, latency)
%   model = fit_latency_model(..., Name=Value)
%
% Inputs, one element per measurement run:
%   location    location label (numbers, strings, cell array of char or categorical)
%   throughput  achieved throughput reported by the receiver (Mbps)
%   latency     mean one-way latency of the run (ms)
%
% Name-value arguments (defaults are the settings of Section 4.3.2):
%   Reference    location with gamma = 0 (default: the first location in sorted order)
%   L0Bounds     [0 50] ms          BetaBounds   [0 5] ms/Mbps
%   ThetaBounds  [0 5000] ms*Mbps   CMax         5000 Mbps (lower bound: max(throughput) + 1)
%   GammaBounds  [-50 50] ms
%
% Output struct:
%   L0, beta, theta, C       estimates
%   offsets                  table of the locations and their offsets gamma_q (ms)
%   MAE, RMSE, R2, n         in-sample fit metrics and number of runs
%   atBound                  names of the parameters that ended at a bound
%   predict(b, location)     predicted latency (ms); predict(b) uses the mean of the
%                            offsets, as for an unmeasured location in the article
%
% The fit is bounded nonlinear least squares (lsqcurvefit, Optimization
% Toolbox) started from several initial values of C; the best fit is kept.
%
% Example (reproduces Table 5 of the article), run from the code folder:
%   T = readtable(fullfile('..', 'data', 'latency_with_achieved_bw_210.csv'));
%   model = fit_latency_model(T.qIdx, T.achieved_Mbps_server, T.latency_ms);
%   model.predict(100, 1)

arguments
    location (:,1)
    throughput (:,1) double {mustBeReal, mustBeFinite}
    latency (:,1) double {mustBeReal, mustBeFinite}
    opts.Reference = []
    opts.L0Bounds (1,2) double = [0 50]
    opts.BetaBounds (1,2) double = [0 5]
    opts.ThetaBounds (1,2) double = [0 5000]
    opts.CMax (1,1) double = 5000
    opts.GammaBounds (1,2) double = [-50 50]
end

n = numel(latency);
if numel(location) ~= n || numel(throughput) ~= n
    error('fit_latency_model:size', 'location, throughput and latency must have the same length.');
end
if iscellstr(location) || iscategorical(location) || ischar(location)
    location = string(location);
end

% location index q = 1..Q, reference location first
locs = unique(location);
if ~isempty(opts.Reference)
    r = find(string(locs) == string(opts.Reference), 1);
    if isempty(r)
        error('fit_latency_model:reference', 'Reference location not found in the data.');
    end
    locs = locs([r, setdiff(1:numel(locs), r)]);
end
[~, q] = ismember(string(location), string(locs));
Q = numel(locs);

bMax = max(throughput);
if bMax + 1 >= opts.CMax
    error('fit_latency_model:CMax', 'CMax must be larger than max(throughput) + 1.');
end
if n <= 4 + (Q-1)
    error('fit_latency_model:data', 'At least %d runs are needed for %d locations.', 4 + Q, Q);
end

lb = [opts.L0Bounds(1), opts.BetaBounds(1), opts.ThetaBounds(1), bMax + 1, ...
      opts.GammaBounds(1)*ones(1, Q-1)];
ub = [opts.L0Bounds(2), opts.BetaBounds(2), opts.ThetaBounds(2), opts.CMax, ...
      opts.GammaBounds(2)*ones(1, Q-1)];
fun = @(p, X) p(1) + p(2).*X(:,1) + p(3)./(p(4) - X(:,1)) + offsetOf(p, X(:,2));
fitOpts = optimoptions('lsqcurvefit', 'Display','off', 'MaxFunctionEvaluations', 2e5, ...
    'MaxIterations', 2e4, 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);

% initial values: for each start value of C the model is linear in the other
% parameters, which are then taken from a linear least-squares fit
C0list = unique(min(opts.CMax, max(bMax + 1, bMax*[1 1.02 1.1 1.25 1.5 2 4])));
D = double(q == 2:Q);                        % location indicators for q = 2..Q
best = inf;
for C0 = C0list
    lin = [ones(n,1), throughput, 1./(C0 - throughput), D] \ latency;
    p0 = [lin(1:3)', C0, lin(4:end)'];
    p0 = min(max(p0, lb), ub);
    [p, resnorm] = lsqcurvefit(fun, p0, [throughput q], latency, lb, ub, fitOpts);
    if resnorm < best
        best = resnorm; pHat = p;
    end
end

gamma = [0, pHat(5:end)]';
Lhat = fun(pHat, [throughput q]);
res = latency - Lhat;

model.L0 = pHat(1);
model.beta = pHat(2);
model.theta = pHat(3);
model.C = pHat(4);
model.offsets = table(locs, gamma, 'VariableNames', {'location', 'gamma_ms'});
model.MAE = mean(abs(res));
model.RMSE = sqrt(mean(res.^2));
model.R2 = 1 - sum(res.^2) / sum((latency - mean(latency)).^2);
model.n = n;
names = [{'L0', 'beta', 'theta', 'C'}, compose('gamma_%s', string(locs(2:end))')];
tol = 1e-6 * max(1, abs(ub - lb));
model.atBound = names(abs(pHat - lb) < tol | abs(ub - pHat) < tol);
model.predict = @(b, varargin) predictLatency(pHat, locs, gamma, b, varargin{:});

if ~isempty(model.atBound)
    warning('fit_latency_model:atBound', ...
        'Parameters at a bound: %s. Consider wider bounds (name-value arguments).', ...
        strjoin(model.atBound, ', '));
end
end

function v = offsetOf(p, q)
g = [0, p(5:end)];
v = g(q)';
v = v(:);
end

function L = predictLatency(p, locs, gamma, b, location)
if any(b >= p(4))
    error('fit_latency_model:capacity', 'Throughput must be below C = %.3f Mbps.', p(4));
end
if nargin < 5
    g = mean(gamma);
else
    [found, k] = ismember(string(location), string(locs));
    if ~all(found)
        error('fit_latency_model:location', 'Unknown location; call predict(b) for an unmeasured location.');
    end
    g = gamma(k);
    g = reshape(g, size(location));
end
L = p(1) + p(2).*b + p(3)./(p(4) - b) + g;
end
