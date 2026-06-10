% Reproduces Table III (PIS, ET, CE under rarer events)
% Paper: Ghazal, Ben Rached, Al-Naffouri (2025), §VI.A
% Setting: M=8, m=4, mu=0.5, gamma_th in {0.4, 0.3, 0.2}
% Runtime: ~1 min on a standard laptop (update after timing)
% Output: table matching Table III in the paper, printed to console
%
% This script handles PIS in both Case 1 (|mu| <= 1) and Case 2 (|mu| > 1)
% of Algorithm 1, so it can also produce Table V (high-mean) rows by
% changing mu and gammaRange at the top.

%% Input

M          = 8;
m          = 4;
mu         = 0.5;
gammaRange = [0.4, 0.3, 0.2];
N_MC       = 1e6;

nG = numel(gammaRange);

% Storage
P_PIS_v    = zeros(1, nG);
RE_PIS_v   = zeros(1, nG);
WNRV_PIS_v = zeros(1, nG);
P_ET_v     = zeros(1, nG);
RE_ET_v    = zeros(1, nG);
WNRV_ET_v  = zeros(1, nG);
P_CE_v     = zeros(1, nG);
RE_CE_v    = zeros(1, nG);
WNRV_CE_v  = zeros(1, nG);

%% PIS

%% PIS (general m, both |mu| <= 1 and |mu| > 1)

q = floor(M / m);

for ig = 1:nG
    gamma_th = gammaRange(ig);

    rng('default')

    NormConstant = ncx2cdf(2*gamma_th, 2*m, 2*m*abs(mu)^2);

    if abs(mu) <= 1
        % Case 1: mode at zero, simple bound
        Ml = gamma_th^m * exp(-m*abs(mu)^2) / (factorial(m) * NormConstant);
        Mg = 1;
    elseif 2*gamma_th < 2*abs(mu)^2 - 2
        % Case 2a: |mu| > 1 and gamma_th below the mode region
        Kmu = ncx2pdf(2*gamma_th, 2, 2*abs(mu)^2);
        Ml  = (2*gamma_th*Kmu)^m / (factorial(m) * NormConstant);
        Mg  = exp(m*abs(mu)^2) * (2*Kmu)^m;
    else
        % Case 2b: |mu| > 1 and gamma_th in/above the mode region
        Lambda = 2*abs(mu)^2;
        x0     = Lambda - 2 + 3/(2*Lambda);
        C1     = exp((Lambda - 1 - x0) * (-0.5 + (sqrt(1 + 4*Lambda*x0) - 1)/(4*x0)));
        C2     = exp((-Lambda + 2 + x0) * (0.5  - (sqrt(1 + Lambda*x0) - 1)/(2*x0)));
        C      = max(C1, C2);
        Kmu    = C * ncx2pdf(x0, 2, 2*abs(mu)^2);
        Ml     = (2*gamma_th*Kmu)^m / (factorial(m) * NormConstant);
        Mg     = exp(m*abs(mu)^2) * (2*Kmu)^m;
    end

    chunk_size      = min(1e7, ceil(1.02 * Ml * q * N_MC));
    unif_chunk_size = chunk_size;

    tic

    P4 = 0;
    i      = 1;
    i_unif = 1;
    buffer      = exprnd(1, m+1, chunk_size);
    unif_buffer = unifrnd(0, 1, unif_chunk_size, 1);
    yM = zeros(m, q);

    for n = 1:N_MC
        for count = 1:q
            while true
                if i > chunk_size
                    buffer = exprnd(1, m+1, chunk_size);
                    i      = 1;
                end
                if i_unif > unif_chunk_size
                    unif_buffer = unifrnd(0, 1, unif_chunk_size, 1);
                    i_unif      = 1;
                end
                Xi     = buffer(:, i);
                Ui     = Xi(1:m) / sum(Xi);
                U      = unif_buffer(i_unif);
                i_unif = i_unif + 1;
                i      = i + 1;
                bound  = prod(exp(-Ui*gamma_th)) / Mg;
                if U < bound || U < prod(besseli(0, 2*abs(mu) .* sqrt(gamma_th*Ui))) * bound
                    break;
                end
            end
            yM(:, count) = Ui;
        end
        y  = reshape(yM, [], 1);
        xo = maxk(y, m);
        if sum(xo) < 1
            P4 = P4 + 1;
        end
    end

    ell_1_PIS = ncx2cdf(2*gamma_th, 2*m, 2*m*mu^2)^q;
    P_PIS     = P4 / N_MC * ell_1_PIS;
    Time_PIS  = toc;

    Var_PIS  = ell_1_PIS * P_PIS - P_PIS^2;
    SCV_PIS  = Var_PIS / P_PIS^2;
    RE_PIS   = sqrt(SCV_PIS / N_MC);
    WNRV_PIS = RE_PIS^2 * Time_PIS;

    P_PIS_v(ig)    = P_PIS;
    RE_PIS_v(ig)   = RE_PIS;
    WNRV_PIS_v(ig) = WNRV_PIS;
end

%% ET

for ig = 1:nG
    gamma_th = gammaRange(ig);

    rng('default')
    tic

    fac  = gamma_th / (2*M);
    S    = 0;
    Term = zeros(N_MC, 1);
    Xb   = fac * (randn(M, N_MC).^2 + randn(M, N_MC).^2);

    for k = 1:N_MC
        X      = Xb(:, k);
        sumOrX = sum(maxk(X, m));
        if sumOrX < gamma_th
            sumX    = sum(X);
            Prod    = prod(besseli(0, 2*abs(mu) * X.^0.5));
            Term(k) = exp((M - gamma_th) * sumX / gamma_th) * Prod;
            S       = S + Term(k);
        end
    end

    P_ET    = S * gamma_th^M * exp(-M*abs(mu)^2) / (M^M * N_MC);
    Term    = Term * gamma_th^M * exp(-M*abs(mu)^2) / (M^M);
    Time_ET = toc;

    Var_ET  = var(Term);
    SCV_ET  = Var_ET / P_ET^2;
    RE_ET   = sqrt(SCV_ET / N_MC);
    WNRV_ET = RE_ET^2 * Time_ET;

    P_ET_v(ig)    = P_ET;
    RE_ET_v(ig)   = RE_ET;
    WNRV_ET_v(ig) = WNRV_ET;
end

%% CE

for ig = 1:nG
    gamma_th = gammaRange(ig);

    rng('default')
    tic

    N_MC_inner = 1e3;
    N_MC_final = N_MC;
    rho        = 0.1;
    low        = floor(rho * N_MC_inner);

    vt      = [0.5, 2*mu^2];
    X       = randn(2*M, N_MC_inner);
    y       = 0.5 * ((X(1:M, :) + mu).^2 + (X(M+1:2*M, :) + mu).^2);
    S       = sum(maxk(y, m));
    [S_so, I] = sort(S);
    S_low   = S_so(1:low);
    gamma_t = S_low(end);
    X_low   = y(:, I(1:low));

    fun = @(v) -(2*vt(1)*exp(vt(2)/2 - mu^2))^M * ...
        sum(exp((0.5/vt(1) - 1)*sum(X_low, 1)) ...
        .* prod(besseli(0, 2*mu*sqrt(X_low)) ./ besseli(0, sqrt(vt(2)/vt(1)*X_low)), 1) ...
        .* (-M*log(2*v(1)) - M*v(2)/2 - sum(X_low, 1)*0.5/v(1) ...
        + sum(log(besseli(0, sqrt(v(2)/v(1)*X_low))), 1)));
    vt      = fminsearch(fun, vt);
    stepsN  = 1;

    while gamma_t > gamma_th
        X       = randn(2*M, N_MC_inner);
        mut     = sqrt(vt(2)/2);
        y       = vt(1) * ((X(1:M, :) + mut).^2 + (X(M+1:2*M, :) + mut).^2);
        S       = sum(maxk(y, m));
        [S_so, I] = sort(S);
        S_low   = S_so(1:low);
        gamma_t = S_low(end);
        if gamma_t < gamma_th
            gamma_t = gamma_th;
            low     = find(S_so > gamma_th, 1, 'first');
        end
        X_low = y(:, I(1:low));
        fun = @(v) -(2*vt(1)*exp(vt(2)/2 - mu^2))^M * ...
            sum(exp((0.5/vt(1) - 1)*sum(X_low, 1)) ...
            .* prod(besseli(0, 2*mu*sqrt(X_low)) ./ besseli(0, sqrt(vt(2)/vt(1)*X_low)), 1) ...
            .* (-M*log(2*v(1)) - M*v(2)/2 - sum(X_low, 1)*0.5/v(1) ...
            + sum(log(besseli(0, sqrt(v(2)/v(1)*X_low))), 1)));
        vt     = fminsearch(fun, vt);
        stepsN = stepsN + 1;
    end

    mut       = sqrt(vt(2)/2);
    X         = randn(2*M, N_MC_final);
    y         = vt(1) * ((X(1:M, :) + mut).^2 + (X(M+1:2*M, :) + mut).^2);
    S         = sum(maxk(y, m));
    [S_so, I] = sort(S);
    low       = find(S_so > gamma_th, 1, 'first');
    X_low     = y(:, I(1:low));

    P_CE   = (2*vt(1)*exp(vt(2)/2 - mu^2))^M * ...
              sum(exp((0.5/vt(1) - 1)*sum(X_low, 1)) ...
              .* prod(besseli(0, 2*mu*sqrt(X_low)) ./ besseli(0, sqrt(vt(2)/vt(1)*X_low)), 1)) / N_MC_final;
    Var_CE = (2*vt(1)*exp(vt(2)/2 - mu^2))^(2*M) * ...
              sum((exp((0.5/vt(1) - 1)*sum(X_low, 1)) ...
              .* prod(besseli(0, 2*mu*sqrt(X_low)) ./ besseli(0, sqrt(vt(2)/vt(1)*X_low)), 1)).^2) / N_MC_final ...
              - P_CE^2;
    Time_CE = toc;

    SCV_CE  = Var_CE / P_CE^2;
    RE_CE   = sqrt(SCV_CE / N_MC_final);
    WNRV_CE = RE_CE^2 * Time_CE;

    P_CE_v(ig)    = P_CE;
    RE_CE_v(ig)   = RE_CE;
    WNRV_CE_v(ig) = WNRV_CE;
end

%% Print Table III

fprintf('\n');
fprintf('TABLE III\n');
fprintf('VARIANCE SUMMARY WITH RARER EVENTS, M = %d, m = %d, mu = %g, N = 10^6\n', M, m, mu);
fprintf('\n');
fprintf('              %-32s    %-32s    %-32s\n', 'PIS', 'ET', 'CE');
fprintf('  %-9s | %-10s  %-8s  %-8s | %-10s  %-8s  %-8s | %-10s  %-8s  %-8s\n', ...
    'gamma_th', 'P_hat_S', 'RE (%)', 'WNRV', 'P_hat_S', 'RE (%)', 'WNRV', 'P_hat_S', 'RE (%)', 'WNRV');
fprintf(['  ----------+--------------------------------+-' ...
         '-------------------------------+--------------------------------\n']);
for ig = 1:nG
    fprintf(['  %-9.4g | %-10.3e  %-8.2f  %-8.1e | %-10.3e  %-8.2f  %-8.1e | ' ...
             '%-10.3e  %-8.2f  %-8.1e\n'], ...
        gammaRange(ig), ...
        P_PIS_v(ig), 100*RE_PIS_v(ig), WNRV_PIS_v(ig), ...
        P_ET_v(ig),  100*RE_ET_v(ig),  WNRV_ET_v(ig), ...
        P_CE_v(ig),  100*RE_CE_v(ig),  WNRV_CE_v(ig));
end
fprintf('\n');