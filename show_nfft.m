% Simple test of the nonequispaced FFT
% 
% Author  Stefan Kunis

% Polynomial degree
N=256;

kk=((-N/2):(N/2-1))';

% Fourier coefficients
f_hat=sin(10*kk/N).*exp(-(kk/N^0.8).^2);

% Showing the Fourier coefficients
subplot(2,2,1);
stem(kk,f_hat);
title('Fourier coefficients');

% Equispaced sampling nodes
xx=linspace(0,1-1/N,N)';

% Clustered sampling nodes
alpha=1/100;
x=xx.^alpha;

% Showing the equispaced sampled values of the trigonometric polynomial
subplot(2,2,2);
semilogy(xx,abs(fft(f_hat)));
ax=axis;
title('Log plot of the FFT');

% Computing nonequispaced DFT and FFT
ft_opt.method='direct';
f1=ndft(f_hat,x,N,ft_opt,'notransp');

ft_opt.method='gaussian';
ft_opt.m=6;                             % N>2*m
ft_opt.sigma=2;
f2=nfft(f_hat,x,N,ft_opt,'notransp');

% Showing the nonequispaced sampled values of the trigonometric polynomial
subplot(2,2,3);
semilogy(xx,abs(fft(f_hat)),x,abs(f1),'o',x,abs(f2),'.');
axis(ax);
title('Log plot of the NFFT');

% Zooming the cluster
subplot(2,2,4);
loglog(xx,abs(fft(f_hat)),'o',x,abs(f1),x,abs(f2),'.');
axis([1-6/N,1,ax(3),ax(4)]);
title('Additionally, logarithm of the nodes');
