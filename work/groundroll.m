% generate dispersive ground roll in the frequency domain
% Use Yedlin (1981) paper
% Daniel Trad - December 2002
echo off
close all
clear d m y M mm;

N=256;
dt=0.04;

m=zeros(N,N);
t=0:2*N-1;
p=-N/2:N/2-1;p=p*2e-4;
t=t*dt;
h=-N/2:N/2-1;h=h*25;

% Once calculated comment out for speed
[F3]=smearingfilt(t,h,p,'LRT',N/4);

M=fft(m,2*N);

% Generate vector with indexes for slopes
maxy=125;
i=1:maxy;
  y(i)=30*(1-exp(-(0.02*(i-1)).^2));

figure(1);
  plot(fix(y),'o');

Nslopes=5;

for ll=1:Nslopes  
  for ii=2:maxy
    jj=10+N-80+ll+round(y(ii-1));
    ii,jj;
    M(ii,jj)=1;
    M(ii+40,jj)=1;
    M(ii+2*40,jj)=1;
    %M(ii,N-jj)=1;
  end
end

window=hanning(maxy+2*40);  
for ii=1:N
  M(1:maxy+2*40,ii)=M(1:maxy+2*40,ii).*window;
end

MM=duplic(M(1:N,:));
figure(2);imagesc(M);colorbar;

% shift the signal down (linear phase shift).
omega=2*pi*freqaxis(dt,length(t));
H=exp(-sqrt(-1)*(omega(:)*ones(1,N)));
mm=ifft(MM.*H);


m=mm(1:N,:);


m=[m;zeros(size(m))];

[d]=applyrtop(m,t,h,F3);
d=real(d);



%a=eye(N);
%mask=cumsum(a);
%mask=[mask;ones(N,N)];
%d=d.*mask;

dd=-1*d(:,N/2+1:N);

dd=[dd zeros(2*N,1)];

figure(3);wigb(dd,2);

writesudata('/home/dtrad/work/groll.bin',dd);







