% Optimized reproduction of the stochastic rows of Table II.
%
% Preprint setting:
%   M = 8, m = 4, mu = 0.5, gamma_th in {1, 0.5}, h_i ~ CN(mu, 1).
%
% Table II sample counts:
%   NMC: N = 1e9 for gamma_th = 1, N = 1e10 for gamma_th = 0.5.
%   MLS: s = 1.6e5 samples per level, 50 independent replicates.
%
% This file keeps the two expensive generic baselines optimized:
%   1) NMC is chunked and fully vectorized within each chunk.
%   2) MLS follows Ben Rached et al. (2020), but carries only the
%      surviving Gamma-process states and vectorizes each level.
% It also appends the remaining stochastic methods used in Table II:
%   UIS, PIS, ET, CE, ETC, and CEC.
%
% Defaults are smoke-test sized so the script can be run quickly. Set
% run_nmc_full and run_mls_full to true to reproduce the paper settings.
% Full time run total execution time: ~ 6 hr on a 

clear; clc;

%% Configuration

M = 8;
m = 4;
mu = 0.5;
gammaRange = [1, 0.5];

run_nmc_full = true;
run_mls_full = true;
run_other_methods = true;

N_NMC_full  = [1e9, 1e10];
N_NMC_smoke = [2e5, 2e5];
N_NMC       = N_NMC_smoke;
if run_nmc_full
    N_NMC = N_NMC_full;
end

N_MLS_full       = 1.6e5;
N_MLS_smoke      = 5e3;
m_split_full     = 50;
m_split_smoke    = 5;
N_MLS            = N_MLS_smoke;
m_split          = m_split_smoke;
if run_mls_full
    N_MLS   = N_MLS_full;
    m_split = m_split_full;
end

N_UIS = 5e6;
N_IS  = 1e6;  % PIS, ET, CE, ETC, CEC

% Larger chunks are usually faster if memory allows. For M = 8, a chunk of
% 2e6 stores only moderate temporary arrays on typical desktops.
nmc_chunk_size = 2e6;

% MLS level-selection heuristic from the 2020 MLS paper. These match the
% attached reproduction script.
s_pilot = 1e3;
pbar    = 0.1;
delta0  = 0.1;

progress_n_updates = 20;
ts = @() char(datetime('now', 'Format', 'HH:mm:ss'));
tic_total = tic;

fprintf('\n[%s] Optimized Table II NMC/MLS run\n', ts());
fprintf('M=%d, m=%d, mu=%g, gammas=%s\n', M, m, mu, mat2str(gammaRange));
fprintf('NMC N=%s, MLS s=%.3g, MLS replicates=%d\n', mat2str(N_NMC), N_MLS, m_split);
fprintf('Other methods: UIS N=%.1e, PIS/ET/CE/ETC/CEC N=%.1e\n', N_UIS, N_IS);
if ~run_nmc_full || ~run_mls_full
    fprintf('Smoke mode is enabled. Set run_nmc_full/run_mls_full true for Table II counts.\n');
end

methods = {'NMC'; 'UIS'; 'PIS'; 'MLS'; 'ET'; 'CE'; 'ETC'; 'CEC'};
nMethods = numel(methods);
P_v     = nan(nMethods, numel(gammaRange));
RE_v    = nan(nMethods, numel(gammaRange));
WNRV_v  = nan(nMethods, numel(gammaRange));
Time_v  = nan(nMethods, numel(gammaRange));
N_v     = nan(nMethods, numel(gammaRange));

Sigma  = eye(M);
SigmaI = eye(M);
rho    = 0;
muV    = mu * ones(M, 1);

%% NMC

fprintf('\n[%s] === Optimized NMC ===\n', ts());
for ig = 1:numel(gammaRange)
    gamma_th = gammaRange(ig);
    rng('default');

    out = run_nmc_optimized(M, m, mu, gamma_th, N_NMC(ig), nmc_chunk_size, ...
        progress_n_updates, ts);

    P_v(1, ig)    = out.P;
    RE_v(1, ig)   = out.RE;
    WNRV_v(1, ig) = out.WNRV;
    Time_v(1, ig) = out.time;
    N_v(1, ig)    = out.N;

    full_seconds = out.seconds_per_sample * N_NMC_full(ig);
    fprintf('[%s] NMC gamma=%g done: P=%.4e, RE=%.2f%%, time=%.1fs.\n', ...
        ts(), gamma_th, out.P, 100*out.RE, out.time);
    fprintf('         observed %.3g samples/s; projected full Table II N=%.1e time: %s.\n', ...
        1/out.seconds_per_sample, N_NMC_full(ig), seconds_to_hms(full_seconds));
end

%% UIS, PIS

if run_other_methods
    fprintf('\n[%s] === UIS ===\n', ts());
    for ig = 1:numel(gammaRange)
        gamma_th = gammaRange(ig);
        rng('default');

        out = run_uis_method(M, m, mu, gamma_th, N_UIS, progress_n_updates, ts);

        P_v(2, ig)    = out.P;
        RE_v(2, ig)   = out.RE;
        WNRV_v(2, ig) = out.WNRV;
        Time_v(2, ig) = out.time;
        N_v(2, ig)    = out.N;
    end

    fprintf('\n[%s] === PIS ===\n', ts());
    for ig = 1:numel(gammaRange)
        gamma_th = gammaRange(ig);
        rng('default');

        out = run_pis_method(M, m, mu, gamma_th, N_IS, progress_n_updates, ts);

        P_v(3, ig)    = out.P;
        RE_v(3, ig)   = out.RE;
        WNRV_v(3, ig) = out.WNRV;
        Time_v(3, ig) = out.time;
        N_v(3, ig)    = out.N;
    end
end

%% MLS

fprintf('\n[%s] === Optimized MLS ===\n', ts());
for ig = 1:numel(gammaRange)
    gamma_th = gammaRange(ig);
    rng('default');

    out = run_mls_optimized(M, m, mu, gamma_th, N_MLS, m_split, ...
        s_pilot, pbar, delta0, progress_n_updates, ts);

    P_v(4, ig)    = out.P;
    RE_v(4, ig)   = out.RE;
    WNRV_v(4, ig) = out.WNRV;
    Time_v(4, ig) = out.mean_replicate_time;
    N_v(4, ig)    = out.N_per_level;

    scale = (N_MLS_full / N_MLS) * (m_split_full / m_split);
    projected_table_time = out.mean_replicate_time * (N_MLS_full / N_MLS);
    projected_total = out.total_primary_time * scale + out.pilot_time;
    fprintf('[%s] MLS gamma=%g done: P=%.4e, RE=%.2f%%, mean replicate time=%.1fs.\n', ...
        ts(), gamma_th, out.P, 100*out.RE, out.mean_replicate_time);
    fprintf('         levels=%s; projected Table-II time/replicate: %s; projected 50-rep wall time: %s.\n', ...
        mat2str(out.t_levels, 4), seconds_to_hms(projected_table_time), seconds_to_hms(projected_total));
end

%% ET, CE, ETC, CEC

if run_other_methods
    fprintf('\n[%s] === ET ===\n', ts());
    for ig = 1:numel(gammaRange)
        gamma_th = gammaRange(ig);
        rng('default');

        out = run_et_method(M, m, mu, gamma_th, N_IS, progress_n_updates, ts);

        P_v(5, ig)    = out.P;
        RE_v(5, ig)   = out.RE;
        WNRV_v(5, ig) = out.WNRV;
        Time_v(5, ig) = out.time;
        N_v(5, ig)    = out.N;
    end

    fprintf('\n[%s] === CE ===\n', ts());
    for ig = 1:numel(gammaRange)
        gamma_th = gammaRange(ig);
        rng('default');

        out = run_ce_method(M, m, mu, gamma_th, N_IS, progress_n_updates, ts);

        P_v(6, ig)    = out.P;
        RE_v(6, ig)   = out.RE;
        WNRV_v(6, ig) = out.WNRV;
        Time_v(6, ig) = out.time;
        N_v(6, ig)    = out.N;
    end

    fprintf('\n[%s] === ETC ===\n', ts());
    for ig = 1:numel(gammaRange)
        gamma_th = gammaRange(ig);
        rng('default');

        out = run_etc_method(M, m, muV, Sigma, SigmaI, gamma_th, N_IS, progress_n_updates, ts);

        P_v(7, ig)    = out.P;
        RE_v(7, ig)   = out.RE;
        WNRV_v(7, ig) = out.WNRV;
        Time_v(7, ig) = out.time;
        N_v(7, ig)    = out.N;
    end

    fprintf('\n[%s] === CEC ===\n', ts());
    for ig = 1:numel(gammaRange)
        gamma_th = gammaRange(ig);
        rng('default');

        out = run_cec_method(M, m, mu, rho, gamma_th, N_IS, progress_n_updates, ts);

        P_v(8, ig)    = out.P;
        RE_v(8, ig)   = out.RE;
        WNRV_v(8, ig) = out.WNRV;
        Time_v(8, ig) = out.time;
        N_v(8, ig)    = out.N;
    end
end

%% Print compact Table II-style rows

fprintf('\nTABLE II rows reproduced by this script\n');
fprintf('M = %d, m = %d, mu = %g\n\n', M, m, mu);

for ig = 1:numel(gammaRange)
    fprintf('gamma_th = %g\n', gammaRange(ig));
    fprintf('  %-8s | %-10s %-12s %-10s %-8s %-10s\n', ...
        'Method', 'N', 'P_hat', 'Time (s)', 'RE (%)', 'WNRV');
    fprintf('  ---------+------------------------------------------------------------\n');
    for im = 1:numel(methods)
        if isnan(P_v(im, ig))
            fprintf('  %-8s | %-10.1e %-12s %-10s %-8s %-10s\n', ...
                methods{im}, N_v(im, ig), '---', '---', '---', '---');
        else
            fprintf('  %-8s | %-10.1e %-12.4e %-10.1f %-8.2f %-10.3e\n', ...
                methods{im}, N_v(im, ig), P_v(im, ig), Time_v(im, ig), ...
                100*RE_v(im, ig), WNRV_v(im, ig));
        end
    end
    fprintf('\n');
end

fprintf('[%s] Finished. Total elapsed time = %s.\n', ts(), seconds_to_hms(toc(tic_total)));


%% Print Table II LaTeX rows

fprintf('\nTABLE II rows reproduced by this script\n');
fprintf('\nTABLE II LaTeX rows\n');
fprintf('M = %d, m = %d, mu = %g\n\n', M, m, mu);

for ig = 1:numel(gammaRange)
    fprintf('gamma_th = %g\n', gammaRange(ig));
    fprintf('  %-8s | %-10s %-12s %-10s %-8s %-10s\n', ...
        'Method', 'N', 'P_hat', 'Time (s)', 'RE (%)', 'WNRV');
    fprintf('  ---------+------------------------------------------------------------\n');
    for im = 1:numel(methods)
        if isnan(P_v(im, ig))
            fprintf('  %-8s | %-10.1e %-12s %-10s %-8s %-10s\n', ...
                methods{im}, N_v(im, ig), '---', '---', '---', '---');
        else
            fprintf('  %-8s | %-10.1e %-12.4e %-10.1f %-8.2f %-10.3e\n', ...
                methods{im}, N_v(im, ig), P_v(im, ig), Time_v(im, ig), ...
                100*RE_v(im, ig), WNRV_v(im, ig));
        end
    end
    fprintf('\n');
[~, ig1] = min(abs(gammaRange - 1));
[~, ig05] = min(abs(gammaRange - 0.5));

for im = 1:numel(methods)
    fprintf('%s & $%s$ &$%s$ &$%.0f$ &$%.1f$ &$%s$    &$%s$ &$%s$ &$%.0f$ &$%.1f$ &$%s$ \\\\ \\hline\n', ...
        methods{im}, ...
        latex_n_value(N_v(im, ig1)), latex_sci_value(P_v(im, ig1), 3), Time_v(im, ig1), 100 * RE_v(im, ig1), latex_sci_value(WNRV_v(im, ig1), 1), ...
        latex_n_value(N_v(im, ig05)), latex_sci_value(P_v(im, ig05), 3), Time_v(im, ig05), 100 * RE_v(im, ig05), latex_sci_value(WNRV_v(im, ig05), 1));
end

fprintf('[%s] Finished. Total elapsed time = %s.\n', ts(), seconds_to_hms(toc(tic_total)));
    out = struct('P', P, 'RE', RE, 'WNRV', WNRV, 'time', elapsed, 'N', N_display);
end

function text = latex_n_value(x)
    if isnan(x)
        text = '---';
        return;
    end

    exponent = floor(log10(x));
    mantissa = x / 10^exponent;
    if abs(mantissa - 1) < 10 * eps(mantissa)
        text = sprintf('10^{%d}', exponent);
    else
        text = sprintf('%.1f\\times 10^{%d}', mantissa, exponent);
    end
end

function text = latex_sci_value(x, decimals)
    if isnan(x)
        text = '---';
        return;
    end
    if x == 0
        text = sprintf('%.*f', decimals, 0);
        return;
    end

    exponent = floor(log10(abs(x)));
    mantissa = x / 10^exponent;
    text = sprintf('%.*f\\times 10^{%d}', decimals, mantissa, exponent);
end

%% Local functions

function out = run_uis_method(M, m, mu, gamma_th, N, progress_n_updates, ts)
    fprintf('\n[%s] UIS gamma=%g: N=%.3g\n', ts(), gamma_th, N);
    timer = tic;

    fff = ncx2cdf(2 * gamma_th, 2, 2 * mu^2);
    chunk_size = min(N, 5e5);
    n_chunks = ceil(N / chunk_size);
    report_stride = max(1, ceil(n_chunks / progress_n_updates));
    hits = 0;

    for c = 1:n_chunks
        this_chunk = min(chunk_size, N - (c - 1) * chunk_size);
        y = 0.5 * ncx2inv(rand(M, this_chunk) * fff, 2, 2 * abs(mu)^2);
        hits = hits + nnz(sum_largest_m(y, m) < gamma_th);

        if c == 1 || mod(c, report_stride) == 0 || c == n_chunks
            done = min(c * chunk_size, N);
            fprintf('[%s]   chunk %d/%d, %.1f%%, conditional hit rate=%.3e, elapsed=%s\n', ...
                ts(), c, n_chunks, 100 * c / n_chunks, hits / done, seconds_to_hms(toc(timer)));
        end
    end

    ell_1 = fff^M;
    P = hits / N * ell_1;
    Var1 = ell_1 * P - P^2;
    out = finalize_mc_out(P, Var1, N, toc(timer), N);
    fprintf('[%s] UIS gamma=%g done: P=%.4e, RE=%.2f%%, time=%.1fs.\n', ...
        ts(), gamma_th, out.P, 100 * out.RE, out.time);
end

function out = run_pis_method(M, m, mu, gamma_th, N, progress_n_updates, ts)
    fprintf('\n[%s] PIS gamma=%g: N=%.3g\n', ts(), gamma_th, N);
    q = floor(M / m);
    norm_constant = ncx2cdf(2 * gamma_th, 2 * m, 2 * m * abs(mu)^2);

    if abs(mu) <= 1
        Ml = gamma_th^m * exp(-m * abs(mu)^2) / (factorial(m) * norm_constant);
        Mg = 1;
    elseif 2 * gamma_th < 2 * abs(mu)^2 - 2
        Kmu = ncx2pdf(2 * gamma_th, 2, 2 * abs(mu)^2);
        Ml  = (2 * gamma_th * Kmu)^m / (factorial(m) * norm_constant);
        Mg  = exp(m * abs(mu)^2) * (2 * Kmu)^m;
    else
        Lambda = 2 * abs(mu)^2;
        x0     = Lambda - 2 + 3 / (2 * Lambda);
        C1     = exp((Lambda - 1 - x0) * (-0.5 + (sqrt(1 + 4 * Lambda * x0) - 1) / (4 * x0)));
        C2     = exp((-Lambda + 2 + x0) * (0.5 - (sqrt(1 + Lambda * x0) - 1) / (2 * x0)));
        Kmu    = max(C1, C2) * ncx2pdf(x0, 2, 2 * abs(mu)^2);
        Ml     = (2 * gamma_th * Kmu)^m / (factorial(m) * norm_constant);
        Mg     = exp(m * abs(mu)^2) * (2 * Kmu)^m;
    end

    chunk_size = min(1e7, max(1, ceil(1.02 * Ml * q * N)));
    fprintf('[%s]   norm_constant=%.3e, Ml=%.3g, Mg=%.3g, buffer=%.3g\n', ...
        ts(), norm_constant, Ml, Mg, chunk_size);

    timer = tic;
    hits = 0;
    i = 1;
    i_unif = 1;
    buffer = -log(rand(m + 1, chunk_size));
    unif_buffer = rand(chunk_size, 1);
    yM = zeros(m, q);
    report_stride = max(1, ceil(N / progress_n_updates));

    for n = 1:N
        for count = 1:q
            while true
                if i > chunk_size
                    buffer = -log(rand(m + 1, chunk_size));
                    i = 1;
                end
                if i_unif > chunk_size
                    unif_buffer = rand(chunk_size, 1);
                    i_unif = 1;
                end

                Xi = buffer(:, i);
                Ui = Xi(1:m) / sum(Xi);
                U = unif_buffer(i_unif);
                i_unif = i_unif + 1;
                i = i + 1;
                bound = prod(exp(-Ui * gamma_th)) / Mg;
                if U < bound || U < prod(besseli(0, 2 * abs(mu) .* sqrt(gamma_th * Ui))) * bound
                    break;
                end
            end
            yM(:, count) = Ui;
        end

        if sum(maxk(reshape(yM, [], 1), m)) < 1
            hits = hits + 1;
        end
        if n == 1 || mod(n, report_stride) == 0 || n == N
            fprintf('[%s]   sample %d/%d, %.1f%%, conditional hit rate=%.3e, elapsed=%s\n', ...
                ts(), n, N, 100 * n / N, hits / n, seconds_to_hms(toc(timer)));
        end
    end

    ell_1 = norm_constant^q;
    P = hits / N * ell_1;
    Var1 = ell_1 * P - P^2;
    out = finalize_mc_out(P, Var1, N, toc(timer), N);
    fprintf('[%s] PIS gamma=%g done: P=%.4e, RE=%.2f%%, time=%.1fs.\n', ...
        ts(), gamma_th, out.P, 100 * out.RE, out.time);
end

function out = run_et_method(M, m, mu, gamma_th, N, progress_n_updates, ts)
    fprintf('\n[%s] ET gamma=%g: N=%.3g\n', ts(), gamma_th, N);
    timer = tic;
    fac = gamma_th / (2 * M);
    Xb = fac * (randn(M, N).^2 + randn(M, N).^2);
    Term = zeros(N, 1);
    report_stride = max(1, ceil(N / progress_n_updates));
    hits = 0;

    for k = 1:N
        X = Xb(:, k);
        if sum(maxk(X, m)) < gamma_th
            hits = hits + 1;
            sumX = sum(X);
            Prod = prod(besseli(0, 2 * abs(mu) * sqrt(X)));
            Term(k) = exp((M - gamma_th) * sumX / gamma_th) * Prod;
        end
        if k == 1 || mod(k, report_stride) == 0 || k == N
            fprintf('[%s]   sample %d/%d, %.1f%%, proposal hit rate=%.3e, elapsed=%s\n', ...
                ts(), k, N, 100 * k / N, hits / k, seconds_to_hms(toc(timer)));
        end
    end

    scale = gamma_th^M * exp(-M * abs(mu)^2) / M^M;
    samples = Term * scale;
    P = mean(samples);
    out = finalize_mc_out(P, var(samples), N, toc(timer), N);
    fprintf('[%s] ET gamma=%g done: P=%.4e, RE=%.2f%%, time=%.1fs.\n', ...
        ts(), gamma_th, out.P, 100 * out.RE, out.time);
end

function out = run_ce_method(M, m, mu, gamma_th, N, progress_n_updates, ts)
    fprintf('\n[%s] CE gamma=%g: adaptive phase, then N=%.3g\n', ts(), gamma_th, N);
    timer = tic;
    N_inner = 1e3;
    rho_CE = 0.1;
    low = floor(rho_CE * N_inner);

    vt = [0.5, 2 * mu^2];
    X = randn(2 * M, N_inner);
    y = 0.5 * ((X(1:M, :) + mu).^2 + (X(M + 1:2 * M, :) + mu).^2);
    S = sum(maxk(y, m));
    [S_so, I] = sort(S);
    gamma_t = S_so(low);
    X_low = y(:, I(1:low));
    steps = 0;

    while true
        steps = steps + 1;
        vt_prev = vt;
        fun = @(v) -(2 * vt_prev(1) * exp(vt_prev(2) / 2 - mu^2))^M * ...
            sum(exp((0.5 / vt_prev(1) - 1) * sum(X_low, 1)) ...
            .* prod(besseli(0, 2 * mu * sqrt(X_low)) ./ besseli(0, sqrt(vt_prev(2) / vt_prev(1) * X_low)), 1) ...
            .* (-M * log(2 * v(1)) - M * v(2) / 2 - sum(X_low, 1) * 0.5 / v(1) ...
            + sum(log(besseli(0, sqrt(v(2) / v(1) * X_low))), 1)));
        vt = fminsearch(fun, vt);
        fprintf('[%s]   CE step %d: vt=[%.4g %.4g], gamma_t=%.4g\n', ts(), steps, vt(1), vt(2), gamma_t);

        if gamma_t <= gamma_th
            break;
        end

        mut = sqrt(vt(2) / 2);
        X = randn(2 * M, N_inner);
        y = vt(1) * ((X(1:M, :) + mut).^2 + (X(M + 1:2 * M, :) + mut).^2);
        S = sum(maxk(y, m));
        [S_so, I] = sort(S);
        gamma_t = S_so(low);
        if gamma_t < gamma_th
            gamma_t = gamma_th;
            low = find(S_so > gamma_th, 1, 'first');
            if isempty(low)
                low = N_inner;
            end
        end
        X_low = y(:, I(1:low));
    end

    mut = sqrt(vt(2) / 2);
    X = randn(2 * M, N);
    y = vt(1) * ((X(1:M, :) + mut).^2 + (X(M + 1:2 * M, :) + mut).^2);
    S = sum(maxk(y, m));
    [S_so, I] = sort(S);
    low = find(S_so > gamma_th, 1, 'first');
    if isempty(low)
        low = N;
    end
    X_low = y(:, I(1:low));
    weights = (2 * vt(1) * exp(vt(2) / 2 - mu^2))^M * ...
        exp((0.5 / vt(1) - 1) * sum(X_low, 1)) ...
        .* prod(besseli(0, 2 * mu * sqrt(X_low)) ./ besseli(0, sqrt(vt(2) / vt(1) * X_low)), 1);
    samples = zeros(N, 1);
    samples(I(1:low)) = weights;
    out = finalize_mc_out(mean(samples), var(samples), N, toc(timer), N);
    fprintf('[%s] CE gamma=%g done: P=%.4e, RE=%.2f%%, time=%.1fs.\n', ...
        ts(), gamma_th, out.P, 100 * out.RE, out.time);
end

function out = run_etc_method(M, m, muV, Sigma, SigmaI, gamma_th, N, progress_n_updates, ts)
    fprintf('\n[%s] ETC gamma=%g: N=%.3g\n', ts(), gamma_th, N);
    timer = tic;
    fac = gamma_th / M;
    facs = sqrt(fac) / sqrt(2);
    Xb = facs * (randn(M, N) + 1i * randn(M, N));
    Term = zeros(N, 1);
    report_stride = max(1, ceil(N / progress_n_updates));
    hits = 0;

    for k = 1:N
        X = Xb(:, k);
        if sum(maxk(abs(X).^2, m)) < gamma_th
            hits = hits + 1;
            Term(k) = exp(real(X' * X / fac - (X - muV)' * SigmaI * (X - muV)));
        end
        if k == 1 || mod(k, report_stride) == 0 || k == N
            fprintf('[%s]   sample %d/%d, %.1f%%, proposal hit rate=%.3e, elapsed=%s\n', ...
                ts(), k, N, 100 * k / N, hits / k, seconds_to_hms(toc(timer)));
        end
    end

    samples = Term * fac^M / det(Sigma);
    out = finalize_mc_out(mean(samples), var(samples), N, toc(timer), N);
    fprintf('[%s] ETC gamma=%g done: P=%.4e, RE=%.2f%%, time=%.1fs.\n', ...
        ts(), gamma_th, out.P, 100 * out.RE, out.time);
end

function out = run_cec_method(M, m, mu, rho, gamma_th, N, progress_n_updates, ts)
    fprintf('\n[%s] CEC gamma=%g: adaptive phase, then N=%.3g\n', ts(), gamma_th, N);
    timer = tic;
    N_inner = 1e3;
    rhoCE = 0.1;
    OM = ones(M, 1);
    JM = ones(M);
    low = floor(rhoCE * N_inner);

    vt = [mu, 1, rho];
    X = (randn(M, N_inner) + 1i * randn(M, N_inner)) / sqrt(2) + mu;
    y = abs(X).^2;
    S = sum(maxk(y, m));
    [S_so, I] = sort(S);
    gamma_t = S_so(low);
    X_low = X(:, I(1:low));
    SigmaRho = vt(2) * ((1 - vt(3)) * eye(M) + vt(3) * JM);
    steps = 0;

    while true
        steps = steps + 1;
        vt_prev = vt;
        fun = @(v) sum( ...
            (M * log(v(2)) + (M - 1) * log(1 - v(3)) + log(1 + (M - 1) * v(3)) ...
            + 1 / (v(2) * (1 - v(3))) * real(sum(conj(X_low - v(1) * OM) .* ...
            ((eye(M) - v(3) / (1 + (M - 1) * v(3)) * JM) * (X_low - v(1) * OM))))) ...
            .* cec_weight_terms(X_low, M, mu, rho, vt_prev, OM, JM));
        vt = fminsearch(fun, vt);
        SigmaRho = vt(2) * ((1 - vt(3)) * eye(M) + vt(3) * JM);
        fprintf('[%s]   CEC step %d: vt=[%.4g %.4g %.4g], gamma_t=%.4g\n', ...
            ts(), steps, vt(1), vt(2), vt(3), gamma_t);

        if gamma_t <= gamma_th
            break;
        end

        X = SigmaRho^0.5 * (randn(M, N_inner) + 1i * randn(M, N_inner)) / sqrt(2) + vt(1);
        y = abs(X).^2;
        S = sum(maxk(y, m));
        [S_so, I] = sort(S);
        gamma_t = S_so(low);
        if gamma_t < gamma_th
            gamma_t = gamma_th;
            low = find(S_so > gamma_th, 1, 'first');
            if isempty(low)
                low = N_inner;
            end
        end
        X_low = X(:, I(1:low));
    end

    X = SigmaRho^0.5 * (randn(M, N) + 1i * randn(M, N)) / sqrt(2) + vt(1);
    y = abs(X).^2;
    S = sum(maxk(y, m));
    [S_so, I] = sort(S);
    low = find(S_so > gamma_th, 1, 'first');
    if isempty(low)
        low = N;
    end
    X_low = X(:, I(1:low));
    weights = cec_weight_terms(X_low, M, mu, rho, vt, OM, JM);
    samples = zeros(N, 1);
    samples(I(1:low)) = weights;
    out = finalize_mc_out(mean(samples), var(samples), N, toc(timer), N);
    fprintf('[%s] CEC gamma=%g done: P=%.4e, RE=%.2f%%, time=%.1fs.\n', ...
        ts(), gamma_th, out.P, 100 * out.RE, out.time);
end

function weights = cec_weight_terms(X_low, M, mu, rho, vt, OM, JM)
    R0 = eye(M) - rho / (1 + (M - 1) * rho) * JM;
    Rt = eye(M) - vt(3) / (1 + (M - 1) * vt(3)) * JM;
    z0 = X_low - mu * OM;
    zt = X_low - vt(1) * OM;

    scale = vt(2)^M * (1 - vt(3))^(M - 1) * (1 + (M - 1) * vt(3)) ...
        / ((1 - rho)^(M - 1) * (1 + (M - 1) * rho));
    exponent = -real(sum(conj(z0) .* (R0 * z0))) / (1 - rho) ...
        + real(sum(conj(zt) .* (Rt * zt))) / (vt(2) * (1 - vt(3)));
    weights = scale * exp(exponent);
end

function out = finalize_mc_out(P, Var1, N_for_re, elapsed, N_display)
    if P > 0 && Var1 >= 0
        RE = sqrt(Var1 / (N_for_re * P^2));
        WNRV = RE^2 * elapsed;
    else
        RE = NaN;
        WNRV = NaN;
    end

    out = struct('P', P, 'RE', RE, 'WNRV', WNRV, 'time', elapsed, 'N', N_display);
end

function out = run_nmc_optimized(M, m, mu, gamma_th, N, chunk_size, progress_n_updates, ts)
    fprintf('\n[%s] NMC gamma=%g: N=%.3g, chunk=%.3g\n', ts(), gamma_th, N, chunk_size);
    timer = tic;

    n_chunks = ceil(N / chunk_size);
    report_stride = max(1, ceil(n_chunks / progress_n_updates));
    mu_shift = sqrt(2) * mu;
    hits = 0;

    for c = 1:n_chunks
        this_chunk = min(chunk_size, N - (c - 1) * chunk_size);

        % If Z ~ CN(mu,1), then |Z|^2 = 0.5*((N(0,1)+sqrt(2)mu)^2+N(0,1)^2).
        y = 0.5 * ((randn(M, this_chunk) + mu_shift).^2 + randn(M, this_chunk).^2);
        sum_top = sum_largest_m(y, m);
        hits = hits + nnz(sum_top < gamma_th);

        if c == 1 || mod(c, report_stride) == 0 || c == n_chunks
            done = min(c * chunk_size, N);
            fprintf('[%s]   chunk %d/%d, %.1f%%, P_hat=%.3e, elapsed=%s\n', ...
                ts(), c, n_chunks, 100 * c / n_chunks, hits / done, seconds_to_hms(toc(timer)));
        end
    end

    elapsed = toc(timer);
    P = hits / N;
    Var1 = P * (1 - P);
    if P > 0
        RE = sqrt(Var1 / (N * P^2));
        WNRV = RE^2 * elapsed;
    else
        RE = NaN;
        WNRV = NaN;
    end

    out = struct('P', P, 'RE', RE, 'WNRV', WNRV, 'time', elapsed, ...
        'N', N, 'seconds_per_sample', elapsed / N, 'hits', hits);
end

function out = run_mls_optimized(M, m, mu, gamma_th, s, m_split, s_pilot, pbar, delta0, progress_n_updates, ts)
    fprintf('\n[%s] MLS gamma=%g: selecting levels with s_pilot=%.3g\n', ts(), gamma_th, s_pilot);
    pilot_timer = tic;
    [t_levels, pilot_estimates] = select_mls_levels(M, m, mu, gamma_th, s_pilot, pbar, delta0, ts);
    pilot_time = toc(pilot_timer);
    L = numel(t_levels);

    fprintf('[%s] MLS gamma=%g: pilot done in %s; t_levels=%s; pilot estimates=%s\n', ...
        ts(), gamma_th, seconds_to_hms(pilot_time), mat2str(t_levels, 4), mat2str(pilot_estimates, 4));

    estimates = zeros(m_split, 1);
    rep_times = zeros(m_split, 1);
    report_stride = max(1, ceil(m_split / progress_n_updates));

    for r = 1:m_split
        rep_timer = tic;
        counts = zeros(1, L);

        G = gamrnd(t_levels(1), 1, M, s);
        keep = outage_from_gamma(G, m, mu, gamma_th);
        counts(1) = nnz(keep);
        survivors = G(:, keep);

        for ell = 2:L
            previous_count = counts(ell - 1);
            if previous_count == 0
                break;
            end

            idx = randi(previous_count, 1, s);
            delta = t_levels(ell) - t_levels(ell - 1);
            G = survivors(:, idx) + gamrnd(delta, 1, M, s);
            keep = outage_from_gamma(G, m, mu, gamma_th);
            counts(ell) = nnz(keep);
            survivors = G(:, keep);
        end

        estimates(r) = prod(counts / s);
        rep_times(r) = toc(rep_timer);

        if r == 1 || mod(r, report_stride) == 0 || r == m_split
            fprintf('[%s]   replicate %d/%d: P=%.3e, counts=%s, time=%s\n', ...
                ts(), r, m_split, estimates(r), mat2str(counts), seconds_to_hms(rep_times(r)));
        end
    end

    P = mean(estimates);
    Var = var(estimates);
    if P > 0
        SCV = Var / P^2;
        RE = sqrt(SCV / m_split);
        WNRV = RE^2 * mean(rep_times);
    else
        RE = NaN;
        WNRV = NaN;
    end

    out = struct('P', P, 'RE', RE, 'WNRV', WNRV, ...
        'mean_replicate_time', mean(rep_times), 'total_primary_time', sum(rep_times), ...
        'pilot_time', pilot_time, 'N_per_level', s, 'm_split', m_split, ...
        't_levels', t_levels, 'replicate_estimates', estimates);
end

function [t_levels, pilot_estimates] = select_mls_levels(M, m, mu, gamma_th, s, pbar, delta0, ts)
    t_grid = [];
    est_grid = [];
    t_prev = 0;
    delta = delta0;
    phat = 1;
    current_states = zeros(M, s);
    level = 1;

    while t_prev < 1
        delta = min(delta, 1 - t_prev);
        accepted = false;

        while ~accepted
            G_try = current_states + gamrnd(delta, 1, M, s);
            keep = outage_from_gamma(G_try, m, mu, gamma_th);
            count = nnz(keep);
            cond_prob = count / s;

            if cond_prob >= pbar || delta <= eps(1)
                accepted = true;
                t_prev = t_prev + delta;
                phat = phat * cond_prob;
                t_grid(end + 1) = t_prev; %#ok<AGROW>
                est_grid(end + 1) = phat; %#ok<AGROW>

                fprintf('[%s]   pilot level %d: t=%.4f, hits=%d/%d, cond=%.3e, est=%.3e\n', ...
                    ts(), level, t_prev, count, s, cond_prob, phat);

                if count == 0 || t_prev >= 1
                    break;
                end

                survivors = G_try(:, keep);
                current_states = survivors(:, randi(count, 1, s));
                delta = min(2 * delta, 1 - t_prev);
                level = level + 1;
            else
                delta = delta / 2;
            end
        end

        if count == 0 || t_prev >= 1
            break;
        end
    end

    est_vec = [1, est_grid];
    time_vec = [0, t_grid];
    min_power = ceil(log10(est_vec(end)));
    phat_targets = 10 .^ (-1:-1:min_power);

    % interp1 is happier when the x-grid is increasing.
    [est_unique, ia] = unique(fliplr(est_vec), 'stable');
    time_unique = fliplr(time_vec);
    time_unique = time_unique(ia);
    t_interp = interp1(est_unique, time_unique, fliplr(phat_targets), 'linear', 'extrap');
    t_levels = [fliplr(t_interp), 1];
    t_levels = unique(max(0, min(1, t_levels)), 'stable');
    pilot_estimates = est_grid;
end

function keep = outage_from_gamma(G, m, mu, gamma_th)
    U = 1 - exp(-G);
    X = 0.5 * ncx2inv(U, 2, 2 * mu^2);
    keep = sum_largest_m(X, m) <= gamma_th;
end

function s = sum_largest_m(X, m)
    if m == size(X, 1)
        s = sum(X, 1);
    else
        s = sum(maxk(X, m, 1), 1);
    end
end

function text = seconds_to_hms(seconds)
    if ~isfinite(seconds)
        text = 'n/a';
        return;
    end

    hours = floor(seconds / 3600);
    minutes = floor((seconds - 3600 * hours) / 60);
    secs = seconds - 3600 * hours - 60 * minutes;

    if hours > 0
        text = sprintf('%dh %02dm %04.1fs', hours, minutes, secs);
    elseif minutes > 0
        text = sprintf('%dm %04.1fs', minutes, secs);
    else
        text = sprintf('%.1fs', secs);
    end
end
