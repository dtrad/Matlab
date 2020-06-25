function [y]=f_k_filter(x)
X=fft2(x);
Y=mute_matrix(X);
y=ifft2(Y);