dt=0.004;
f=35;
N=256;
ispike=200;
x=zeros(N,1);
x(ispike,1)=1;
w=ricker(f,dt);
xx=conv(x,w);x=xx(1:N);
t=0:N-1;t=t*dt;

figure(1)
[S,freq,xi]=stransform0(x,t);
subplot(311);imagesc(abs(S));colorbar;
subplot(312);plot(x);
subplot(313);plot(xi);

y=nmo(x,dt,1000,1500);
df=strfilter(x,dt,1000,1500);

if (0) 
  figure(2);
  [S,freq,yi]=stransform0(y,t);
  subplot(311);imagesc(abs(S));colorbar;
  subplot(312);plot(y);
  subplot(313);plot(yi);
end

figure(3);
[S,freq,yi]=stransform_filter(y,t,df);
subplot(311);imagesc(abs(S));colorbar;
subplot(312);plot(y);
subplot(313);plot(yi);