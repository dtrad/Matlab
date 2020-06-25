M=512;
x=randn(1,M);
x=[x 0 x(M:-1:2)];
X=toeplitz(x);
y=randn(1,M);
y=[y zeros(size(y))];

figure

subplot(211);plot(X*y.');
subplot(212);plot(real(ifft(fft(x).*fft(y))))
