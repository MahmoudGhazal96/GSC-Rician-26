% Reproduces Table VII (ETC + CEC, constant correlation)
% Paper: Ghazal, Ben Rached, Al-Naffouri (2025), §VI.B
% Runtime: ~20 sec on a standard laptop
% Output: table matching Table VII in the paper, printed to console

%% Input

M    = 8;
m    = 4;
mu   = 0.5;
rho  = 0.6;
gammaRange = [0.075, 0.1, 0.2, 0.3, 0.4];

Sigma = (1-rho)*eye(M) + rho*ones(M);
muV   = mu*ones(M,1);

% Storage for results
nG       = numel(gammaRange);
P_ET_v   = zeros(1,nG);
RE_ET_v  = zeros(1,nG);
WNRV_ET_v = zeros(1,nG);
P_CE_v   = zeros(1,nG);
RE_CE_v  = zeros(1,nG);
WNRV_CE_v = zeros(1,nG);

for ig = 1:nG
    gamma = gammaRange(ig);

    %% ET

    rng('default')
    N_MC = 1e6;
    tic
    fac    = gamma/M;
    facs   = fac^0.5/sqrt(2);
    S      = 0;
    Term   = zeros(N_MC,1);
    Xb     = facs*(randn(M,N_MC) + 1i*randn(M,N_MC));
    SigmaI = Sigma^-1;
    for k = 1:N_MC
        X      = Xb(:,k);
        Xn     = abs(Xb(:,k)).^2;
        sumOrX = sum(maxk(Xn,m));
        if sumOrX < gamma
            Term(k) = exp(X'*X/fac - (X-muV)'*SigmaI*(X-muV));
            S       = S + Term(k);
        end
    end

    P_ET    = S*fac^M/(det(Sigma)*N_MC);
    Term    = Term*fac^M/det(Sigma);
    Time_ET = toc;

    ell      = P_ET;
    Var_ET   = var(Term);
    SCV_ET   = Var_ET/ell^2;
    RE_ET    = (SCV_ET/N_MC)^0.5;
    WNRV_ET  = RE_ET^2*Time_ET;

    P_ET_v(ig)    = P_ET;
    RE_ET_v(ig)   = RE_ET;
    WNRV_ET_v(ig) = WNRV_ET;

    %% CE

    rng('default')
    tic
    N_MC    = 1e3;
    N_MC_1  = 1e6;
    rhoCE   = 0.1;
    OM      = ones(M,1);
    JM      = ones(M);
    low     = floor(rhoCE*N_MC);

    vt = [mu,1,rho];
    X  = (randn(M,N_MC) + 1i*randn(M,N_MC))/sqrt(2) + mu;
    y  = abs(X).^2;
    S  = sum(maxk(y,m));
    [S_so,I] = sort(S);
    S_low    = S_so(1:low);
    gamma_t  = S_low(end);
    X_low    = X(:,I(1:low));
    fun = @(v) sum((M*log(v(2)) + (M-1)*log(1-v(3)) + log(1+(M-1)*v(3)) ...
        + 1/(v(2)*(1-v(3)))*real(sum(conj(X_low-v(1)*OM).*((eye(M)-v(3)/(1+(M-1)*v(3))*JM)*(X_low-v(1)*OM))))) ...
        .*(vt(2)^M*(1-vt(3))^(M-1)*(1+(M-1)*vt(3))/((1-rho)^(M-1)*(1+(M-1)*rho)) ...
        *exp(-1/(1-rho) ...
        *real(sum(conj(X_low-mu*OM).*((eye(M)-rho/(1+(M-1)*rho)*JM)*(X_low-mu*OM)))) ...
        + 1/(vt(2)*(1-vt(3)))* ...
        real(sum(conj(X_low-vt(1)*OM).*((eye(M)-vt(3)/(1+(M-1)*vt(3))*JM)*(X_low-vt(1)*OM))))) ));
    vt       = fminsearch(fun,vt);
    SigmaRho = vt(2)*((1-vt(3))*eye(M) + vt(3)*JM);
    stepsN   = 1;

    while gamma_t > gamma
        X = SigmaRho^0.5*(randn(M,N_MC) + 1i*randn(M,N_MC))/sqrt(2) + vt(1);
        y = abs(X).^2;
        S = sum(maxk(y,m));
        [S_so,I] = sort(S);
        S_low    = S_so(1:low);
        gamma_t  = S_low(end);
        if gamma_t < gamma
            gamma_t = gamma;
            low     = find(S_so > gamma, 1, 'first');
        end
        X_low = X(:,I(1:low));
        fun = @(v) sum((M*log(v(2)) + (M-1)*log(1-v(3)) + log(1+(M-1)*v(3)) ...
            + 1/(v(2)*(1-v(3)))*real(sum(conj(X_low-v(1)*OM).*((eye(M)-v(3)/(1+(M-1)*v(3))*JM)*(X_low-v(1)*OM))))) ...
            .*(vt(2)^M*(1-vt(3))^(M-1)*(1+(M-1)*vt(3))/((1-rho)^(M-1)*(1+(M-1)*rho)) ...
            *exp(-1/(1-rho) ...
            *real(sum(conj(X_low-mu*OM).*((eye(M)-rho/(1+(M-1)*rho)*JM)*(X_low-mu*OM)))) ...
            + 1/(vt(2)*(1-vt(3)))* ...
            real(sum(conj(X_low-vt(1)*OM).*((eye(M)-vt(3)/(1+(M-1)*vt(3))*JM)*(X_low-vt(1)*OM))))) ));
        vt       = fminsearch(fun,vt);
        SigmaRho = vt(2)*((1-vt(3))*eye(M) + vt(3)*JM);
        stepsN   = stepsN + 1;
    end

    X        = SigmaRho^0.5*(randn(M,N_MC_1) + 1i*randn(M,N_MC_1))/sqrt(2) + vt(1);
    y        = abs(X).^2;
    S        = sum(maxk(y,m));
    [S_so,I] = sort(S);
    low      = find(S_so > gamma, 1, 'first');
    S_low    = S_so(1:low);
    X_low    = X(:,I(1:low));

    P_CE = sum(vt(2)^M*(1-vt(3))^(M-1)*(1+(M-1)*vt(3))/((1-rho)^(M-1)*(1+(M-1)*rho)) ...
            *exp(-1/(1-rho) ...
            *real(sum(conj(X_low-mu*OM).*((eye(M)-rho/(1+(M-1)*rho)*JM)*(X_low-mu*OM)))) ...
            + 1/(vt(2)*(1-vt(3)))* ...
            real(sum(conj(X_low-vt(1)*OM).*((eye(M)-vt(3)/(1+(M-1)*vt(3))*JM)*(X_low-vt(1)*OM))))))/N_MC_1;
    Var_CE = sum(vt(2)^(2*M)*(1-vt(3))^(2*M-2)*(1+(M-1)*vt(3))^2/((1-rho)^(2*M-2)*(1+(M-1)*rho)^2) ...
            *exp(-2/(1-rho) ...
            *real(sum(conj(X_low-mu*OM).*((eye(M)-rho/(1+(M-1)*rho)*JM)*(X_low-mu*OM)))) ...
            + 2/(vt(2)*(1-vt(3)))* ...
            real(sum(conj(X_low-vt(1)*OM).*((eye(M)-vt(3)/(1+(M-1)*vt(3))*JM)*(X_low-vt(1)*OM))))))/N_MC_1 - P_CE^2;

    Time_CE  = toc;

    RE_CE    = (Var_CE/(P_CE^2*N_MC_1))^0.5;
    WNRV_CE  = RE_CE^2*Time_CE;

    P_CE_v(ig)    = P_CE;
    RE_CE_v(ig)   = RE_CE;
    WNRV_CE_v(ig) = WNRV_CE;
end

%% Print Table VII

fprintf('\n');
fprintf('TABLE VII\n');
fprintf('CEC VS ETC UNDER CONSTANT CORRELATION, M = %d, m = %d, mu = %g, rho = %g, N = 10^6\n', M, m, mu, rho);
fprintf('\n');
fprintf('              %-32s    %-32s\n', 'ETC', 'CEC');
fprintf('  %-9s | %-10s  %-8s  %-8s | %-10s  %-8s  %-8s\n', ...
    'gamma_th', 'P_hat_S', 'RE (%)', 'WNRV', 'P_hat_S', 'RE (%)', 'WNRV');
fprintf('  ----------+--------------------------------+---------------------------------\n');
for ig = 1:nG
    fprintf('  %-9.4g | %-10.3e  %-8.2f  %-8.1e | %-10.3e  %-8.2f  %-8.1e\n', ...
        gammaRange(ig), ...
        P_ET_v(ig), 100*RE_ET_v(ig), WNRV_ET_v(ig), ...
        P_CE_v(ig), 100*RE_CE_v(ig), WNRV_CE_v(ig));
end
fprintf('\n');
