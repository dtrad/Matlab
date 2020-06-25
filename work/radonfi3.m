% Program to study the behavior of the smearing filter for the Radon Transform
% The smearing filter is given by one row of LH*L, where L and LH are the 
% projection operators (back and forth) between Radon space and offset space,
% both in frequency. 
% L=exp(-i*w*q*h^2)
% LH=exp(i*w*q*h^2) 
% Such that d=L*m and m_adjoint=LH*d
% Because LH*L is a circulant matrix, LH*L*m is equivalent to a convolution 
% of a row of LH*L with m (see Strang, Introduction to Applied Mathematics) 
% Daniel Trad
% UBC- July 1999

clear
l=-1000:10:1000;
%PP=5*(rand(size(l))-.5);
%l=l+PP;
dq=1e-5;
PMAX=5e-3;
PMIN=-PMAX;
p1=PMIN:dq:PMAX;
p2=0;
w=2*pi*10;
for ii=1:length(p1);
   G(ii)=sum(exp(-i*w*(p2-p1(ii))*l));
end;
GK=abs(fft(G));
kpaxis=1:round(length(p1)/2);
Nyquist=1/(2*w*pi*dq);
kpaxis=kpaxis*Nyquist/max(kpaxis);
subplot(311),plot(p1,real(G)), title('(a)');xlabel('p');ylabel('real(G)');
%subplot(212),plot(kpaxis,GK(1:round(length(p1)/2))); title('(b)');xlabel('kq');ylabel('Amplitude(G)');
clear
l=-1000:10:1000;
PP=0*(rand(size(l))-.5);
l=l+PP;
dq=1e-5;
PMAX=5e-3;
PMIN=-PMAX;
p1=PMIN:dq:PMAX;

p2=0;
w=2*pi*25;
for ii=1:length(p1);
   G(ii)=sum(exp(-i*w*(p2-p1(ii))*l));
end;
GK=abs(fft(G));

kpaxis=1:round(length(p1)/2);
Nyquist=1/(2*w*pi*dq);
kpaxis=kpaxis*Nyquist/max(kpaxis);
subplot(312),plot(p1,real(G)), title('(b)');xlabel('p');ylabel('real(G)');
%subplot(224),plot(kpaxis,GK(1:round(length(p1)/2))); title('(d)');xlabel('kq')%;ylabel('amplitude(G)');


clear
l=-1000:10:1000;
PP=10*(rand(size(l))-.5);
l=l+PP;
dq=1e-5;
PMAX=5e-3;
PMIN=-PMAX;
p1=PMIN:dq:PMAX;

p2=0;
w=2*pi*25;
for ii=1:length(p1);
   G(ii)=sum(exp(-i*w*(p2-p1(ii))*l));
end;
GK=abs(fft(G));

kpaxis=1:round(length(p1)/2);
Nyquist=1/(2*w*pi*dq);
kpaxis=kpaxis*Nyquist/max(kpaxis);
subplot(313),plot(p1,real(G)), title('(c)');xlabel('p');ylabel('real(G)');
%subplot(224),plot(kpaxis,GK(1:round(length(p1)/2))); title('(d)');xlabel('kq')%;ylabel('amplitude(G)');








