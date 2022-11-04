
w=ricker(0.004,2,2)
figure(1);
subplot(121),plot(w)
ww=conv(w,w(end:-1:1),'same');
subplot(122),plot(ww)
figure(2)
t=1:length(w);
plot(t,w,t,ww,'.-')
