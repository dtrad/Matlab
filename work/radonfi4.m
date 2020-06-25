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

l=-2000:40:2000; % Offset axis

dq=0.5e-8; % radon interval
PMAX=10e-7;  
PMIN=-PMAX;
p1=PMIN:dq:PMAX;
p2=0;
w=2*pi*5; % Frequency
for ii=1:length(p1);
   G(ii)=sum(exp(-i*w*(p2-p1(ii))*l.^2));
end;
GK=abs(fft(G));
kpaxis=1:round(length(p1)/2);
Nyquist=1/(2*w*pi*dq);
kpaxis=kpaxis*Nyquist/max(kpaxis);
subplot(311),plot(p1,real(G)), title('(a)');xlabel('q');ylabel('real(G)');
%subplot(222),plot(kpaxis,GK(1:round(length(p1)/2))); title('(b)');xlabel('kq');ylabel('Amplitude');
 
clear

l=-2000:40:2000; % Offset axis

dq=0.5e-8; % radon interval
PMAX=10e-7;  
PMIN=-PMAX;
p1=PMIN:dq:PMAX;
p2=0;
w=2*pi*10; % Frequency
for ii=1:length(p1);
   G(ii)=sum(exp(-i*w*(p2-p1(ii))*l.^2));
end;
GK=abs(fft(G));
kpaxis=1:round(length(p1)/2);
Nyquist=1/(2*w*pi*dq);
kpaxis=kpaxis*Nyquist/max(kpaxis);
subplot(312),plot(p1,real(G)), title('(b)');xlabel('q');ylabel('real(G)');
%subplot(222),plot(kpaxis,GK(1:round(length(p1)/2))); title('(b)');xlabel('kq');ylabel('Amplitude');

l=-2000:40:2000;
PP=40*(rand(size(l))-.5);
l=l+PP;
dq=0.5e-8;
PMAX=10e-7;
PMIN=-PMAX;
p1=PMIN:dq:PMAX;

p2=0;
%w=2*pi*25;
for ii=1:length(p1);
   G(ii)=sum(exp(-i*w*(p2-p1(ii))*l.^2));
end;
GK=abs(fft(G));
kpaxis=1:round(length(p1)/2);
Nyquist=1/(2*w*pi*dq);
kpaxis=kpaxis*Nyquist/max(kpaxis);
subplot(313),plot(p1,real(G)), title('(c)');xlabel('q');ylabel('real(G)');

%subplot(224),plot(kpaxis,GK(1:round(length(p1)/2))); title('(d)');xlabel('kq');ylabel('Amplitude');









