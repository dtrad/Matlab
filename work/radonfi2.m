clear
l=-1000:50:1000;
dq=1e-5;
PMAX=1e-3;
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
subplot(221),plot(p1,real(G)), title('(a)');xlabel('q');ylabel('real(G)');
subplot(222),plot(kpaxis,GK(1:round(length(p1)/2))); title('(b)');xlabel('kq');ylabel('amplitude(G)');
clear
l=-1000:200:1000;
dq=1e-5;
PMAX=1e-3;
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
subplot(223),plot(p1,real(G)), title('(c)');xlabel('q');ylabel('real(G)');
subplot(224),plot(kpaxis,GK(1:round(length(p1)/2))); title('(d)');xlabel('kq');ylabel('amplitude(G)');









