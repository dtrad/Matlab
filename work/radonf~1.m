l=-1000:10:1000;
p1=-1e-6:1e-8:1e-6;
p2=0;
w=2*pi*20;
for ii=1:length(p1);
   G(ii)=sum(exp(-i*w*(p1(ii)-p2)*l.^2));
end;
subplot(211),plot(p1,real(G))
subplot(212),plot(fftshift(abs(fft(G))));
