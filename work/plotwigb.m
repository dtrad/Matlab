function plotwigb
clf
figure(1)

load ridge_result
load wavelet_result
subplot(221);wigb(d0);title('original');

subplot(222);wigb(d);title('noise (order 1)');
subplot(223);wigb(d2);title('wavelet denoising');
subplot(224);wigb(d2r);title('ridgelet denoising');