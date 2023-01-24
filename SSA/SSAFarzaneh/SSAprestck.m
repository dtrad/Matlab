clc, clear all, close all;
load data2

dt=0.001;
k=2; %k:rank of the matrix (Number of events) 
Nitr=5; %Nitr: number of iteration to converge the amplitude
[N,n]=size(datakill);
win=17; %to applty SSA on pre-stack datawe need to assume data is linear by generating small windows
%(for post-stack data win=n)
nwin=floor(n)/(win); lent=win; mid=(win+1)/2;
fS=1/dt; fN=fS/2; fhigh=fS;

fftdata=fft(datakill);
tic
for ii=1:nwin
mm=fftdata(:,(ii-1)*win+1:ii*win);       
for i=1:Nitr   
for j=1:N
i1=(j-1)*fS/(N-1);
d=mm(j,:);

%% function for rank reduction:
[svdh]=rankred(d,k);

%% funtion for extracting the antidiagonals:
[s]=antidiag(svdh,mid,lent);
  
%% interpolation part:
T=d~=0;
iT=ones(size(T));
Sob=d+((iT-T).*s);          
mm(j,:)=Sob;

end
end
fftdata(:,(ii-1)*win+1:ii*win)=mm;
end
toc

newdata=real(ifft(fftdata));
plotseismogram(newdata)
