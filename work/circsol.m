% Solve the deconvolution problem using svd for circulant matrices
% Daniel Trad - April 5 2000

x=rickerm(1,.1);
nt=length(x);
t=1:nt-1;
W=toeplitz(x);
W=W(1:nt-1,1:nt-1);
sf=(reverse(sort(abs(real(fft(W(:,1)))))));
ss=(svd(W));
figure,
semilogy(t,sf,'+',t,ss)


