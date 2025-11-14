
% TOY_SVGD_FWI
% Simple Stein Variational Gradient Descent example mimicking Bayesian FWI.
%
% - m: 1D model (e.g., velocity at n grid points)
% - d: synthetic seismic data from linear forward operator G
% - Posterior: p(m|d) ∝ exp( - 0.5||Gm - d||^2/sigma_d^2 - 0.5||m - m0||^2/sigma_m^2 )
%
% This uses the same setup as the toy HMC example, but applies SVGD
% with an ensemble of particles instead of a Markov chain.
close all
clear
% Set default figure position
screenSize = get(0, 'ScreenSize');
set(0, 'DefaultFigurePosition', [50, screenSize(4)-400-200, 600-100, 400-100]);

rng(0);  % for reproducibility

%% 1. Define "true" model and linear forward operator (toy FWI)
n  = 20;                 % number of model parameters
m_true = 2 + 0.5*exp(-((1:n)' - 10).^2/10);  % smooth bump

% Forward operator G (toy "wave-equation" sensitivity matrix)
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

%% 5. SVGD parameters
Np        = 50;       % number of particles (models)
n_iter    = 2000;     % number of SVGD iterations
eps_step  = 0.02;     % step size for SVGD

% Initialize particles around prior mean
particles = repmat(m0, 1, Np) + 0.3*randn(n, Np);

%% 6. SVGD iterations
for it = 1:n_iter
    % Compute gradients of log posterior at each particle
    grad_log_p = zeros(n, Np);
    for j = 1:Np
        % log p(m) = -U(m) + const, so ∇ log p = -∇U
        grad_log_p(:,j) = -gradU(particles(:,j));
    end

    % Compute RBF kernel matrix and its gradients
    % Pairwise squared distances between particles
    dists2 = zeros(Np, Np);
    for i = 1:Np
        for j = 1:Np
            diff_ij = particles(:,i) - particles(:,j);
            dists2(i,j) = diff_ij' * diff_ij;
        end
    end

    % Bandwidth h using median heuristic
    upper_triu_idx = find(triu(ones(Np),1));
    dist_vec = dists2(upper_triu_idx);
    med_sq = median(dist_vec);
    if med_sq <= 0
        h = 1.0;
    else
        h = med_sq / log(Np + 1);
    end

    K = exp(-dists2 / h);    % RBF kernel matrix

    % SVGD update for each particle
    phi = zeros(n, Np);

    for i = 1:Np
        phi_i = zeros(n,1);
        for j = 1:Np
            % k(x_j, x_i)
            k_ji = K(j,i);

            % Attraction term: k(x_j,x_i) * grad_log_p(x_j)
            attr = k_ji * grad_log_p(:,j);

            % Repulsion term: ∇_{x_j} k(x_j,x_i)
            % For k(x,z) = exp(-||x-z||^2 / h):
            % ∇_x k = -2/h * k * (x - z)
            diff_j = particles(:,j) - particles(:,i);
            rep = (-2.0 / h) * k_ji * diff_j;

            phi_i = phi_i + attr + rep;
        end
        phi(:,i) = phi_i / Np;
    end

    % Update particles
    particles = particles + eps_step * phi;

    % (Optional) print progress occasionally
    if mod(it, 500) == 0
        fprintf('SVGD iteration %d / %d\n', it, n_iter);
        %plot(particles);figure(gcf);pause
    end
end

%% 7. Posterior statistics from particles
m_mean = mean(particles, 2);
m_std  = std(particles, 0, 2);

%% 8. (Optional) Compute MAP solution (Gaussian posterior => closed form)
A = (G'*G)/(sigma_d^2) + (1/sigma_m^2)*eye(n);
b = (G'*d_obs)/(sigma_d^2) + (1/sigma_m^2)*m0;
m_map = A \ b;

%% 9. Plot results
figure; clf;
xgrid = 1:n;

% Left: models
subplot(1,2,1); hold on; box on;
% plot a few particles
plot(xgrid, particles(:,1:10), 'Color', [0.7 0.7 0.7]);
plot(xgrid, m_true, 'k-', 'LineWidth', 2);
plot(xgrid, m_map,  'b--', 'LineWidth', 1.5);
plot(xgrid, m_mean, 'r-', 'LineWidth', 1.5);
legend('Samples (subset)', 'True model', 'MAP (L2 FWI)', ...
       'Posterior mean (SVGD)', 'Location','Best');
xlabel('Model index'); ylabel('m');
title('Toy FWI: True vs MAP vs SVGD posterior mean');

% Right: uncertainty
subplot(1,2,2); hold on; box on;
plot(xgrid, m_std, 'r-', 'LineWidth', 1.5);
xlabel('Model index'); ylabel('Std(m)');
title('Posterior standard deviation (SVGD)');


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
