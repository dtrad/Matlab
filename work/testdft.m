t=0:16;
dt=0.01;
t=t*dt;
t2=-length(t):length(t);
t2=t2*dt;
f1=0;
f2=10;
x=1+sin(2*pi*f1*t)+sin(2*pi*f2*t);
h1=sinc(2*pi*f2*(t2-dt/2));
h2=sinc(2*pi*f2*(t2));

y1=conv(x,h1);
y2=conv(x,h2);
y=[y1(1),y(2)];
for ii=2:length(y1);
  y=[y,y1(ii),y2(ii)];
end  

subplot(311);plot(t,x,'o');
subplot(312);plot(t2,h1);
subplot(313);plot(y,'o');
figure(gcf);
