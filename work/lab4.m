% Lab # 4

nt=128;
f1=20;
f2=60;
dt=0.005;
t=0:nt-1;
t=t*dt;
x=sin(2*pi*f1*t)+sin(2*pi*f2*t);

f=-nt/2:nt/2-1;
f=f*2/nt*1/(2*dt);

figure(1)
subplot(311),plot(t,x);title('x')
subplot(312),plot(f,fftshift(abs(fft(x))));
subplot(313),plot(f,fftshift(angle(fft(x))));

h=[1 0.5];
hz=[h zeros(1,126)];

figure(2)
subplot(211),plot(f,fftshift(abs(fft(hz))));title('h')
subplot(212),plot(f,fftshift(angle(fft(hz))));

y=conv(x,h);y=y(1:nt);

figure(3)
subplot(311),plot(t,y);title('y=conv(h,x)')
subplot(312),plot(f,fftshift(abs(fft(y))));
subplot(313),plot(f,fftshift(angle(fft(y))));

% Part 2
%a)
y=ifft(fft(x).*fft(hz));
figure(4)
subplot(311),plot(t,real(y));title('y conv in freq')
subplot(312),plot(f,fftshift(abs(fft(y))));
subplot(313),plot(f,fftshift(angle(fft(y))));


%b)
h=[0.5 1];
hz=[h zeros(1,126)];

figure(5)
subplot(211),plot(f,fftshift(abs(fft(hz))));title('non minimum phase')
subplot(212),plot(f,fftshift(angle(fft(hz))));

y=conv(x,h);y=y(1:nt);

figure(6)
subplot(311),plot(t,y);title('y=conv(x,h2)')
subplot(312),plot(f,fftshift(abs(fft(y))));
subplot(313),plot(f,fftshift(angle(fft(y))));





