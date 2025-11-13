% ising_binary_inversion_demo.m
% 1-D Ising-style binary inversion with simulated annealing (no toolboxes)

rng(7);

%% 1) True binary model (±1)
N = 64;
m_true = -ones(N,1);
m_true(20:40) = 1;   % "plume-like" block

%% 2) Ricker wavelet and convolution matrix W
dt = 0.004;
w  = ricker(25.0, dt, 17);    % zero-phase, short band-limited
L  = numel(w);
half = floor(L/2);

W = zeros(N,N);
for i = 1:N
    for k = 1:L
        j = i - half + (k-1);
        if j>=1 && j<=N
            W(i,j) = w(k);
        end
    end
end

%% 3) Data with noise
d_clean = W*m_true;
noise = 0.05 * max(abs(d_clean)) * randn(N,1);
d = d_clean + noise;

%% 4) Energy (data misfit + Ising pairwise coupling)
lam = 1.0;                     % pairwise weight
WTW = W.'*W;
WTd = W.'*d;

energy = @(m) (m.'*(WTW*m) - 2*(WTd.'*m)) + ...
              lam * sum(1 - m(1:end-1).*m(2:end));

%% 5) Simulated annealing (random spin flips)
T  = 3.0; T_min = 1e-3; alpha = 0.97; itersPerT = 400;
m0 = sign(rand(N,1)-0.5); m0(m0==0) = 1;

m  = m0;
E  = energy(m);
accepted = 0;

while T > T_min
    for it = 1:itersPerT
        i = randi(N);
        m_new = m; m_new(i) = -m_new(i);     % flip one spin
        E_new = energy(m_new);
        dE = E_new - E;
        if dE <= 0 || rand < exp(-dE/T)
            m = m_new; E = E_new;
            accepted = accepted + 1;
        end
    end
    T = T*alpha;
end

%% 6) Results & plots
fprintf('Initial energy: %.2f\n', energy(m0));
fprintf('Final energy:   %.2f\n', energy(m));
fprintf('Accepted flips: %d\n', accepted);

figure('Position',[100 100 900 350]);
plot(m_true,'LineWidth',1.8); hold on;
plot(m,'--','LineWidth',1.8);
xlabel('Index'); ylabel('Spin (±1)');
title('Ising-style binary inversion (1D)');
legend('True binary model','Estimated (anneal)'); grid on;

d_est = W*m;
figure('Position',[100 100 900 350]);
plot(d,'LineWidth',1.8); hold on;
plot(d_est,'--','LineWidth',1.8);
xlabel('Index'); ylabel('Amplitude');
title('Data fit after inversion');
legend('Observed data','Reconstructed (W m_{est})'); grid on;

%% --------- helper: Ricker wavelet ----------
function w = ricker(f0, dt, nt)
    t = ((0:nt-1) - (nt-1)/2) * dt;
    a = (pi*f0*t).^2;
    w = (1 - 2*a).*exp(-a);
end
