% Program synthrad.m
% Produces a synthetic parabolic radon model and the corresponding data
% Daniel Trad - UBC - EOS - 15-7-99
clear
%close all
curve='HRT';
savesu='y'
pathfile='/home/dtrad/work/';
file='model5';
filename=[pathfile file];

nt=256;
nh=64;
np=40;
dh=25;
dt=0.004;
h_near=-800;
h=0:nh-1;h=h*dh;h=h+h_near;
t=0:nt-1;t=t*dt;
f=50;
Vmin=1200;
Vmax=3000;
theta=45;
ww=ones(nt*np,1);
V(1)=2700;tt(1)=0.2;coef(1)=0.5;
V(2)=3000;tt(2)=0.4;coef(2)=0.5;
V(3)=5000;tt(3)=0.6;coef(3)=-0.5;
%V(4)=1500;tt(4)=0.4;coef(4)=-0.25;
%V(5)=1500;tt(5)=0.6;coef(5)=0.125;
%V(6)=1550;tt(6)=0.2;coef(6)=0.1;


[p,dp]=radonaxis2(np,Vmin,Vmax,h,curve);
%dp=p(2)-p(1);
v=zeros(nt,np);
w=rickerm(f,dt);     
for ii=1:length(V)
   pp(ii)=1/V(ii)^2;
   ip(ii)=round(pp(ii)/dp)+1
   it(ii)=round(tt(ii)/dt)+1
end

for ii=1:length(V) v(it(ii),ip(ii))=coef(ii);end

[nt np]=size(v)
for i=1:np,v2(:,i)=conv(v(:,i),w);end

v=v2(1:nt,1:np);
figure,wigb(v,1,p,t);
if (curve=='PRT')
   d=radonL(v,p,dt,h);
elseif (curve=='HRT')
   if (theta==0) d=radonop(v,t,h,p);
   else d=radonopd(v,t,h,p,theta);
   end
end      
figure,wigb(d,1,h,t);   
if (savesu=='y')
fid = fopen([filename '.bin'],'wb')
for it=1:length(d(:))
         fwrite(fid,d(it),'float32');
end
fclose(fid)
   
fid = fopen([filename '.off'],'wb')
for it=1:length(h)
      fwrite(fid,h(it),'float32');
end
fclose(fid)
end;









