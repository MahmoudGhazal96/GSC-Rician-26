%% =====================================================================
%  Outage probability of GSC(M,L)/MRC in i.i.d. Rician fading
%  via MGF + Laplace inversion (Annamalai et al., IEEE TVT 2006)
%
%  Hyperparameters set as suggested in the paper:
%    Laplace inversion (Abate-Whitt):  A = 30, B = 18, C = 24
%
%  Note: the paper's Table I gives a closed-form marginal MGF in terms of
%  the Marcum-Q function, but MATLAB's built-in marcumq does not accept
%  complex arguments. Since the Laplace inversion evaluates the MGF at
%  complex s, we compute the marginal MGF
%        phi(s,x) = int_x^inf exp(-s t) p(t) dt
%  by numerical quadrature, which handles complex s natively.
% =====================================================================

%% ---- Driver -----------------------------------------------------------
gamma_star = 1;          % threshold SNR
M = 4;                     % number of strongest paths combined
L = 8;                     % total number of diversity paths
mu = 0.5;                  % LOS amplitude

K     = abs(mu)^2;         % Rice factor
Omega = 1 + abs(mu)^2;     % total per-branch power

fprintf('M = %d, m = %d, mu = %.2f, gamma = %.1e\n', L, M, mu, gamma_star);
tic
Pout = gsc_rician_cdf(gamma_star, M, L, K, Omega);
toc
fprintf('Estimated outage probability = %.3e\n', Pout);


%% ---- Outage CDF of GSC(M,L) in Rician fading --------------------------
function Pout = gsc_rician_cdf(gamma_th, M, L, K, Omega)
% Compute outage CDF of GSC(M,L) receiver in i.i.d. Rician fading
% following Annamalai et al. (IEEE TVT 2006).

if ~(isscalar(M) && isscalar(L) && M >= 1 && L >= M)
    error('M,L must be integers with 1 <= M <= L');
end

combML    = nchoosek(L, M);
pdf_fun   = @(t) rician_pdf_snr(t, K, Omega);
phi_marg  = @(s, x) marginal_phi_numeric(s, x, pdf_fun);
phi_gamma = @(s) phi_gamma_numeric(s, M, L, combML, pdf_fun, phi_marg);

% Laplace inversion constants -- paper values
A = 30;
B = 18;
C = 24;

Pout = laplace_inversion_abate_whitt(phi_gamma, gamma_th, A, B, C);
end


%% ---- Rician SNR PDF (per branch) --------------------------------------
function p = rician_pdf_snr(t, K, Omega)
p   = zeros(size(t));
pos = (t >= 0);
if any(pos)
    tp   = t(pos);
    z    = 2 * sqrt(K * (1 + K) .* tp / Omega);
    logp = log(1 + K) - log(Omega) - K - (1 + K) .* tp / Omega ...
           + log(besseli(0, z));
    p(pos) = exp(logp);
end
end


%% ---- Marginal MGF: phi(s,x) = int_x^inf exp(-s t) p(t) dt -------------
function phi = marginal_phi_numeric(s, x, pdf_fun)
% Handles complex s natively. Upper integration limit caps the tail of
% p(t); exp(-(1+K) t) decays fast enough that Tmax = 25 is safe for
% K, Omega of order unity.
Tmax      = 25;
integrand = @(t) exp(-s .* t) .* pdf_fun(t);
opts      = {'RelTol', 1e-6, 'AbsTol', 1e-10, 'ArrayValued', true};
phi       = integral(integrand, x, Tmax, opts{:});
end


%% ---- GSC output MGF (eq. (4) of the paper) ----------------------------
function val = phi_gamma_numeric(s, M, L, combML, pdf_fun, phi_marg)
% Implements eq. (4): integrate over theta in [0, pi/2) via tan substitution.
% Upper limit taken as atan(Tmax) since the integrand decays exponentially
% in t = tan(theta).
Tmax      = 15;
theta_max = atan(Tmax);
integrand = @(theta) inner(theta, s, M, L, combML, pdf_fun, phi_marg);
opts      = {'RelTol', 1e-5, 'AbsTol', 1e-8, 'ArrayValued', true};
val       = integral(integrand, 0, theta_max, opts{:});
end

function y = inner(theta, s, M, L, combML, pdf_fun, phi_marg)
t      = tan(theta);
sec2   = 1 ./ (cos(theta).^2);
n      = numel(t);
pvals  = pdf_fun(t);
phi0   = zeros(1, n);
phi_s  = zeros(1, n);
for k = 1:n
    phi0(k)  = phi_marg(0, t(k));
    phi_s(k) = phi_marg(s, t(k));
end
y = (M * combML) * exp(-s .* t) .* pvals ...
    .* ((1 - phi0).^(L - M)) .* (phi_s.^(M - 1)) .* sec2;
end


%% ---- Laplace inversion (Abate-Whitt; eq. (6) of the paper) ------------
function F = laplace_inversion_abate_whitt(phi_gamma, x, A, B, C)
alpha = @(b) (b == 0) * 0.5 + (b ~= 0);
pref  = 2^(1 - C) * exp(A / 2);
outer = 0;
for c = 0:C
    inner_sum = 0;
    for b = 0:(c + B)
        s         = (A + 1i * 2 * pi * b) / (2 * x);
        denom     = (A + 1i * 2 * pi * b);
        phi       = phi_gamma(s);
        inner_sum = inner_sum + (-1)^b * alpha(b) * (phi / denom);
    end
    outer = outer + nchoosek(C, c) * inner_sum;
end
F = pref * real(outer);
F = max(0, min(1, F));
end