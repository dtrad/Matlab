
% TOY_HMC_FWI
% Simple Hamiltonian Monte Carlo example mimicking Bayesian FWI.
% - m: 1D model (e.g., velocity at n grid points)
% - d: synthetic seismic data from linearized forward operator G
% - Posterior: p(m|d) ∝ exp( - 0.5||Gm - d||^2/sigma_d^2 - 0.5||m - m0||^2/sigma_m^2 )
%
% In a real FWI:
%   G  ~ linearization of wave-equation operator (Jacobian)
%   gradU(m) ~ adjoint-state gradient
% Here we use a small matrix G for simplicity.

rng(0);  % for reproducibility

%% 1. Define "true" model and linear forward operator (toy FWI)
n  = 20;                 % number of model parameters
m_true = 2 + 0.5*exp(-((1:n)' - 10).^2/10);  % smooth bump

% Forward operator G (toy "wave-equation" sensitivity matrix)
% Here: G is a banded smoothing + random matrix to mimic coupling.
Nd = 40;                 % number of data samples
x  = linspace(0,1,n);
t  = linspace(0,1,Nd);
G  = zeros(Nd, n);
for i = 1:Nd
    % Gaussian sensitivity around some position (like kernels)
    G(i,:) = exp(-(x - t(i)).^2 / 0.02);
end
G = G + 0.1*randn(size(G));   % add some variability

%% 2. Generate synthetic noisy data
sigma_d = 0.05;                     % noise std in data
d_clean = G * m_true;
d_obs   = d_clean + sigma_d*randn(Nd,1);

%% 3. Define Gaussian prior on m (e.g., smooth background)
m0       = 2*ones(n,1);   % prior mean
sigma_m  = 0.5;           % prior std

%% 4. Define potential U(m) = -log p(m|d) and its gradient
U = @(m) potential_U(m, G, d_obs, m0, sigma_d, sigma_m);
gradU = @(m) grad_potential_U(m, G, d_obs, m0, sigma_d, sigma_m);

%% 5. HMC parameters
eps_step   = 0.005;   % leapfrog step size
L_steps    = 50;      % number of leapfrog steps per proposal
Nsamples   = 4000;    % total number of HMC samples
burnin     = 1000;    % discard initial samples as burn-in

m_current  = m0;      % start from prior mean

samples    = zeros(n, Nsamples);
accepted   = 0;

%% 6. Run HMC sampling
for k = 1:Nsamples
    % Sample momentum p ~ N(0, I)
    p_current = randn(n,1);

    % Save current state
    m_prop = m_current;
    p_prop = p_current;

    % --- Leapfrog integration ---
    p_prop = p_prop - 0.5*eps_step * gradU(m_prop);
    for i = 1:L_steps
        % position update
        m_prop = m_prop + eps_step * p_prop;

        % momentum update (full step except at last iteration)
        if i < L_steps
            p_prop = p_prop - eps_step * gradU(m_prop);
        end
    end
    p_prop = p_prop - 0.5*eps_step * gradU(m_prop);
    % Negate momentum for symmetry (optional but common)
    p_prop = -p_prop;

    % --- Metropolis acceptance step ---
    H_current = U(m_current) + 0.5*(p_current'*p_current);
    H_prop    = U(m_prop)    + 0.5*(p_prop'*p_prop);
    dH = H_prop - H_current;

    if log(rand) < -dH
        % Accept
        m_current = m_prop;
        accepted  = accepted + 1;
    end

    samples(:,k) = m_current;
end

accept_rate = accepted / Nsamples;
fprintf('HMC acceptance rate: %.2f\n', accept_rate);

%% 7. Discard burn-in, compute posterior mean and std
samples_post = samples(:, burnin+1:end);
m_mean = mean(samples_post, 2);
m_std  = std(samples_post, 0, 2);

%% 8. (Optional) Compute MAP solution (Gaussian posterior => closed form)
A = (G'*G)/(sigma_d^2) + (1/sigma_m^2)*eye(n);
b = (G'*d_obs)/(sigma_d^2) + (1/sigma_m^2)*m0;
m_map = A \ b;

%% 9. Plot results
figure; clf;
xgrid = 1:n;

% Model: true, MAP, posterior mean
subplot(1,2,1);
plot(xgrid, m_true, 'k-', 'LineWidth', 2); hold on;
plot(xgrid, m_map,  'b--', 'LineWidth', 1.5);
plot(xgrid, m_mean, 'r-', 'LineWidth', 1.5);
legend('True model', 'MAP (L2 FWI)', 'Posterior mean (HMC)', ...
       'Location','Best');
xlabel('Model index'); ylabel('m');
title('Toy FWI: True vs MAP vs Posterior mean');
grid on;

% Uncertainty (posterior std)
subplot(1,2,2);
plot(xgrid, m_std, 'r-', 'LineWidth', 1.5);
xlabel('Model index'); ylabel('Std(m)');
title('Posterior standard deviation (uncertainty)');
grid on;


%% --------- Helper functions: potential and gradient -------------------

function Uval = potential_U(m, G, d_obs, m0, sigma_d, sigma_m)
    % Negative log posterior (up to a constant)
    % U(m) = 0.5*||Gm - d||^2 / sigma_d^2 + 0.5*||m - m0||^2 / sigma_m^2
    res_d = G*m - d_obs;
    res_m = m - m0;
    Uval = 0.5*(res_d'*res_d)/(sigma_d^2) + 0.5*(res_m'*res_m)/(sigma_m^2);
end

function g = grad_potential_U(m, G, d_obs, m0, sigma_d, sigma_m)
    % Gradient of U(m) w.r.t. m
    % ∇U = G'*(Gm - d)/sigma_d^2 + (m - m0)/sigma_m^2
    g = (G'*(G*m - d_obs))/(sigma_d^2) + (m - m0)/(sigma_m^2);
end
