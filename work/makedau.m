t = (1:4096)*2;
wave = makewavelet(4,8,'Daubechies',4,'Father',4096);
	subplot(111);plot(t,wave);title(' D4 Wavelet');xlabel('time (s)');
%	plot(t(300:800),wave(300:800)); title(' D4 Wavelet ');
figure(gcf)