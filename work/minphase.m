function [y]=minphase(x)
% Computes the min. phase function from a function x
% the MP function is y
%		[y]=minphase(x)
% Theory: Tad Ulrych: Notes from GEOP520b- UBC-CA
% Daniel Trad- UBC- 30-07-98
x=paddyad(x);
X=fft(x);
B=log(abs(X));
BH=hilbert(B);
Y=exp(BH);
%ly=length(Y);Y=duplic(Y(1:ly/2));
y=ifft(Y);
ly=length(y);
y=y(ly:-1:1);