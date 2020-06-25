clear;
close all;
n1=1024;
n2=128;

load data2m.dat
data2m=data2m.';
[nt nx]=size(data2m);



p0=data2m;

nt2=round(nt/2);
nx2=round(nx/2);

wt2=zeros(n1,n2);
	
wt2(1:nt2,1:nx2)=p0(1:nt2,1:nx2);
wt2(n1-nt2+1:n1,1:nx2)=p0(nt2+1:nt,1:nx2);
wt2(1:nt2,n2-(nx-nx2-1):n2)=p0(1:nt2,nx2+1:nx);
wt2(n1-nt2+1:n1,n2-(nx-nx2-1):n2)=p0(nt2+1:nt,nx2+1:nx);

