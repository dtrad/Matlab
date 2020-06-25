function [X]=delft_dft(x,t);
% dft from delft group;
Xmax=2*pi;

%----- DFT matrix -----
N=length(x);
M=N;

m=(-(M/2):(-1+M/2))';
delta_kx=2*pi/Xmax;

F=exp(i*m*t*delta_kx);

%----- Actual DFT -----
X=F*x(:);
return;