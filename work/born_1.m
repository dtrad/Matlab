% Multiple removal
load Roxt;
dh=10;dt=0.004;
c0=1000;rho=1000;
D=Ro;
D=seis_shape(D);
[nt,nh]=size(D);
w=frequency(dt,nt);
kx=frequency(dh,nh);

D=fft2(D);
kx=frequency(dt,nh);
xs=0;
zs=10;
xg=(1:nh)*dh;
zg=(1:nh)*10;

AA=sin(w(:))*ones(1,length(kx));
BB=((w(:)/c0).^2*ones(1,length(kx)));
CC=ones(length(w),1)*(kx(:).').^2;

q=AA.*(BB-CC).^2;

%*sqrt((w*ones(1,size(kx))/c0).^2-ones(size(w),1)*kx.^2);

G0d=greend(xs,zs,kx,zg,q,rho);
G0fs=greend(xs,-zs,kx,zg,q,rho);
G0=G0d+G0fs;
figure,
G0dt=ifft2(G0d);
G0fst=ifft2(G0fs);
subplot(211);wigb(real(G0dt));
subplot(212);wigb(real(G0fst));


