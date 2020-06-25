% sinc function in Lumley paper 

k=3;
dt=0.004;
dt=0.004;
w=-124:124;
w=2*pi*w;
G=((sin(w*dt*(k+1)/2)).^2)./((sin(w*dt/2)).^2);
figure(1);
plot(G);
figure(gcf);


%filter in t
n = 5;
k = (n-1)/2;
a=[ones(1,n)];

N=256;
t=1:N;
m=(length(a)-1)/2;
x=rand(1,N);
x=sin(2*pi*t*10*dt)+0.3*sin(2*pi*t*120*dt);

y=conv(x,a);
z=y(k:k+N-1);
 

y=conv(z(end:-1:1),a);
z=y(k+N-1:-1:k);

z=z/((2*k)*(2*k));

    

%z=deconv(y,b);
figure(2);
subplot(211);plot(t,x,t,z);
f=freqaxis(dt,N);
subplot(212);plot(f,fftshift(abs(fft(x))),f,fftshift(abs(fft(z))));


% straight convolution
a=ones(1,n);
y=conv(x,a);
zc=y(k+1:k+N);
zc=zc/(n-2);
figure(3)
subplot(211);plot(t,x,t,zc);
subplot(212);plot(f,fftshift(abs(fft(x))),f,fftshift(abs(fft(zc))));

%comparison with Claerbout convbox
yc=boxcarClaerbout(n,N,x);


zc2=boxcarRec(n,N,x);
figure(4)
subplot(111);plot(t,x,t,zc,t,yc,t,zc2,'-o');

return;
figure(3)
k=0;plot(fftshift(abs(fft([-1 zeros(1,k) 2 zeros(1,k) -1 zeros(1,100)]))));pause;
k=1;plot(fftshift(abs(fft([-1 zeros(1,k) 2 zeros(1,k) -1 zeros(1,100)]))));pause;
k=2;plot(fftshift(abs(fft([-1 zeros(1,k) 2 zeros(1,k) -1 zeros(1,100)]))));pause;
k=3;plot(fftshift(abs(fft([-1 zeros(1,k) 2 zeros(1,k) -1 zeros(1,100)]))));pause;


    