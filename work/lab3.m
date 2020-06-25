clear
close all
t=-10:.1:9.9;
tc=-19.9:.1:19.9;
H=[zeros(1,100),ones(1,100)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fa=conv(H,H);
subplot(311),plot(t,H);
subplot(312),plot(t,H);
subplot(313),plot(tc,fa)

x=exp(-t).*H;
fb=conv(H,x);
figure,
subplot(311),plot(t,H)
subplot(312),plot(t,x)
subplot(313),plot(tc,fb)

x=exp(-t).*H;
y=exp(-2*t).*H;
fc=conv(x,y);
figure,
subplot(311),plot(t,x)
subplot(312),plot(t,y)
subplot(313),plot(tc,fc)

%%%%%%%%%%%%%
% e

x=exp(-2*t).*H;
y=exp(-abs(t));
fe=conv(x,y);
figure,
subplot(311),plot(t,x)
subplot(312),plot(t,y)
subplot(313),plot(tc,fe)

%%%%%%%%%%%%%%%%%%%%%%%
delta=zeros(1,100);
delta(50)=1;
subplot(211),plot(abs(fft(delta))),axis([0 100 0 1]);
subplot(212),plot(angle(fft(delta)))
