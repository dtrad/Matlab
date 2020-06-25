function [c]=resta(a,b,plotopt,text)
% Substract vectors of different size;
%  [c]=resta(a,b,plot)
% Daniel Trad - UBC-

if nargin<3 plotopt='n';end;
if nargin<4 text='residuals';end;
la=length(a);;lb=length(b);lc=min(la,lb);

a=a(1:lc);a=a(:);
b=b(1:lc);b=b(:);

c=a-b;
if plotopt=='y'
figure,
subplot(221);plot(c);title(text);
subplot(222);plot(a);title('a');vv=axis;
subplot(223);plot(b);title('b');axis(vv);
subplot(221);axis(vv);
end