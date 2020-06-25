function [z]=circonv(x,y)
% Circular convolution between x and y using fft
% Daniel Trad- UBC- 23-08-98
X=fft(x);
Y=fft(y);
Z=X.*Y;
z=ifft(Z);

