% Reproduces Table VIII (ETC vs. Song's asymptotic, NMC corroboration)
% Paper: Ghazal, Ben Rached, Al-Naffouri (2025), §VI.B
% Distribution: Jakes' correlation, M=8, m=4, K=0.4, Omega=1, d=lambda
% Runtime: ~3 min on a standard laptop (NMC is the dominant cost)
% Output: table matching Table VI in the paper, printed to console
%
% Asymptotic formula from:
%   X. Song and J. Cheng, "Asymptotic Analysis of GSC Over Arbitrarily
%   Correlated Rician Channels," WOCC 2013, eq. (26).
%
% Paper notation -> code notation:
%   N  (total branches)              -> M
%   L  (combined branches)           -> m
%   gamma_th (outage threshold)      -> gamma
%   bar_gamma_i (avg branch SNR)     -> avgSNR  (i.i.d. across branches)
%   K_i (Rician K-factor)            -> K       (i.i.d. across branches)
%   R   (scattering covariance mat.) -> Sigma
%   M   (normalised correlation mat.)-> M_norm  (= Sigma/sigma_sc^2)
%   c_L (LOS component vector)       -> muV     (real, zero phases assumed)
%   A = R^{-1}                       -> SigmaI

%% Input

M         = 8;
m         = 4;
K         = 0.4;
W         = 1;                  % d/lambda (antenna separation / wavelength)
avgSNR    = 10^(10/10);         % linear scale  (= 10 dB)
sigma_sc  = sqrt(avgSNR / (K + 1));   % scattering std
mu        = sqrt(K / (K + 1) * avgSNR);
gammaRange = [1, 2.5, 5, 7.5, 10];
nmc_mask   = [false, false, true, true, true];  % which gammas get NMC

muV  = mu * ones(M, 1);

% Jakes' covariance: R_kl = sigma_sc^2 * J0(2pi|k-l|W)
r     = @(h) (sigma_sc^2) * besselj(0, 2*pi*abs(h)*W);
lags  = repmat((1:M).', 1, M) - repmat(1:M, M, 1);
Sigma = r(lags);
SigmaI = Sigma^-1;
SigmaH = Sigma^0.5;

nG = numel(gammaRange);

% Storage
P_ET_v    = zeros(1, nG);
RE_ET_v   = zeros(1, nG);
WNRV_ET_v = zeros(1, nG);
P_NMC_v   = nan(1, nG);
RE_NMC_v  = nan(1, nG);
WNRV_NMC_v = nan(1, nG);
P_asymp_v = zeros(1, nG);

%% ETC (exponential twisting, correlated)

for ig = 1:nG
    gamma = gammaRange(ig);

    rng('default')
    N_MC = 1e6;
    tic

    fac  = gamma / M;
    facs = sqrt(fac) / sqrt(2);

    S    = 0;
    Term = zeros(N_MC, 1);

    Xb = facs * (randn(M, N_MC) + 1i*randn(M, N_MC));

    for k = 1:N_MC
        X      = Xb(:, k);
        Xn     = abs(X).^2;
        sumOrX = sum(maxk(Xn, m));
        if sumOrX < gamma
            Term(k) = exp(real(X'*X/fac - (X - muV)'*SigmaI*(X - muV)));
            S       = S + Term(k);
        end
    end

    P_ET    = S * fac^M / (det(Sigma) * N_MC);
    Term    = Term * fac^M / det(Sigma);
    Time_ET = toc;

    Var_ET  = var(Term);
    SCV_ET  = Var_ET / P_ET^2;
    RE_ET   = sqrt(SCV_ET / N_MC);
    WNRV_ET = RE_ET^2 * Time_ET;

    P_ET_v(ig)    = P_ET;
    RE_ET_v(ig)   = RE_ET;
    WNRV_ET_v(ig) = WNRV_ET;
end

%% Naive Monte-Carlo (NMC) -- batched, only for the three larger gammas

N_MC  = 1e8;
chunk = 1e6;
if mod(N_MC, chunk) ~= 0
    error('N_MC must be a multiple of chunk.');
end
n_chunks = N_MC / chunk;

for ig = 1:nG
    if ~nmc_mask(ig)
        continue
    end
    gamma = gammaRange(ig);

    rng('default')
    tic

    S = 0;
    for c = 1:n_chunks
        Xb      = SigmaH * (randn(M, chunk) + 1i*randn(M, chunk))/sqrt(2) + muV;
        Xn      = abs(Xb).^2;
        Xn_sort = sort(Xn, 1, 'descend');
        sumTop  = sum(Xn_sort(1:m, :), 1);
        S       = S + sum(sumTop < gamma);
    end
    Time_NMC = toc;

    P_NMC = S / N_MC;
    if P_NMC > 0
        Var_NMC  = P_NMC * (1 - P_NMC);
        SCV_NMC  = Var_NMC / P_NMC^2;
        RE_NMC   = sqrt(SCV_NMC / N_MC);
        WNRV_NMC = RE_NMC^2 * Time_NMC;
    else
        RE_NMC   = NaN;
        WNRV_NMC = NaN;
    end

    P_NMC_v(ig)    = P_NMC;
    RE_NMC_v(ig)   = RE_NMC;
    WNRV_NMC_v(ig) = WNRV_NMC;
end

%% Song's asymptotic approximation -- eq. (26) of [Song & Cheng 2013]

M_norm    = Sigma / sigma_sc^2;
quad_form = real(muV.' * SigmaI * muV);

binom_NL  = nchoosek(M, m);
Gamma_NmL = factorial(M - m);
exp_LOS   = exp(-quad_form);
K_factor  = (1 + K)^M;
fact_N    = factorial(M);
L_pow     = m^(M - m);
det_Mnorm = det(M_norm);
SNR_prod  = avgSNR^M;

for ig = 1:nG
    gth_pow      = gammaRange(ig)^M;
    P_asymp_v(ig) = (binom_NL * Gamma_NmL * exp_LOS * K_factor * gth_pow) / ...
                    (fact_N * L_pow * det_Mnorm * SNR_prod);
end

%% Relative difference (ETC vs Song)

RelDiff_v = 100 * abs(P_asymp_v - P_ET_v) ./ P_ET_v;

%% Print Table VIII

fprintf('\n');
fprintf('TABLE VIII\n');
fprintf(['ETC VS. SONG''S ASYMPTOTIC APPROXIMATION UNDER JAKES'' CORRELATION, ' ...
         'M = %d, m = %d, K = %g, Omega = %g, d = lambda\n'], M, m, K, 1);
fprintf('\n');
fprintf('              %-32s %-12s %-32s    %s\n', 'ETC', 'Song', 'NMC', '');
fprintf('  %-9s | %-10s  %-8s  %-8s | %-10s | %-10s  %-8s  %-8s | %s\n', ...
    'gamma_th', 'P_hat_S', 'RE (%)', 'WNRV', 'P_hat', 'P_hat_S', 'RE (%)', 'WNRV', 'RelDiff (%)');
fprintf(['  ----------+--------------------------------+------------' ...
         '+---------------------------------+------------\n']);
for ig = 1:nG
    if nmc_mask(ig)
        nmc_p  = sprintf('%-10.3e', P_NMC_v(ig));
        nmc_re = sprintf('%-8.2f',  100*RE_NMC_v(ig));
        nmc_w  = sprintf('%-8.1e',  WNRV_NMC_v(ig));
    else
        nmc_p  = sprintf('%-10s', '-');
        nmc_re = sprintf('%-8s',  '-');
        nmc_w  = sprintf('%-8s',  '-');
    end
    fprintf('  %-9.4g | %-10.3e  %-8.2f  %-8.1e | %-10.3e | %s  %s  %s | %.0f\n', ...
        gammaRange(ig), ...
        P_ET_v(ig), 100*RE_ET_v(ig), WNRV_ET_v(ig), ...
        P_asymp_v(ig), ...
        nmc_p, nmc_re, nmc_w, ...
        RelDiff_v(ig));
end
fprintf('\n');
