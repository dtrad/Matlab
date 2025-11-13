% ising_with_ga.m  (requires Global Optimization Toolbox)
rng(7);

N = 64;
m_true = -ones(N,1); m_true(20:40)=1;

dt = 0.004; w = ricker(25.0, dt, 17);
L = numel(w); half = floor(L/2);

W = zeros(N,N);
for i = 1:N
    for k = 1:L
        j = i - half + (k-1);
        if j>=1 && j<=N, W(i,j) = w(k); end
    end
end

d_clean = W*m_true;
noise   = 0.05 * max(abs(d_clean)) * randn(N,1);
d = d_clean + noise;

lam = 1.0; WTW = W.'*W; WTd = W.'*d;

% Objective on x in {0,1}^N, with m = 2x - 1
obj = @(x) energy_ising_cont(2*x(:)-1, WTW, WTd, lam);

opts = optimoptions('ga','PopulationType','bitstring', ...
    'MaxGenerations',200,'Display','iter');

[x_best, fval] = ga(obj, N, [], [], [], [], [], [], [], opts);
m_est = 2*x_best(:)-1;

fprintf('GA objective: %.2f\n', fval);

figure('Position',[100 100 900 350]);
plot(m_true,'LineWidth',1.8); hold on;
plot(m_est,'--','LineWidth',1.8);
xlabel('Index'); ylabel('Spin (±1)');
title('Ising-style binary inversion (GA bitstring)');
legend('True','Estimated'); grid on;

d_est = W*m_est;
figure('Position',[100 100 900 350]);
plot(d,'LineWidth',1.8); hold on;
plot(d_est,'--','LineWidth',1.8);
xlabel('Index'); ylabel('Amplitude');
title('Data fit (GA result)'); legend('Observed','W m_{est}'); grid on;

function val = energy_ising_cont(m, WTW, WTd, lam)
    misfit = m.'*(WTW*m) - 2*(WTd.'*m);
    pair   = sum(1 - m(1:end-1).*m(2:end));
    val = misfit + lam*pair;
end

function w = ricker(f0, dt, nt)
    t = ((0:nt-1) - (nt-1)/2) * dt;
    a = (pi*f0*t).^2;
    w = (1 - 2*a).*exp(-a);
end
