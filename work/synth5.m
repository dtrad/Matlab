% model
clear;
for ii=1:15,close,end;
options='y'; %save
optionm='n'; %multiples
%optionw='minr'; % wavelet
optionw='rick';
optionn='y'; % noise;

tf=1;
dt=0.001;
dh=5;
nh=256;
h_near=0;
nw=50;
f=40;
aten=-1.2;

v(1)=5000;t0(1)=0.05;cr(1)=0.3;
v(2)=5500;t0(2)=0.1;cr(2)=0.3;
v(3)=6000;t0(3)=0.15;cr(3)=-0.2;
v(4)=7000;t0(4)=0.2;cr(4)=0.2;

% multiples
if optionm=='y'
v(5)=v(1);t0(5)=2*t0(1);cr(5)=-cr(1)^2;
v(6)=v(1);t0(6)=3*t0(1);cr(6)=+cr(1)^3;
v(7)=v(1);t0(7)=4*t0(1);cr(7)=-cr(1)^4;
v(8)=v(1);t0(8)=5*t0(1);cr(8)=+cr(1)^5;
v(9)=v(2);t0(9)=t0(3)+t0(1);cr(9)=-cr(3)*cr(1);
v(10)=v(2);t0(10)=t0(3)+2*t0(1);cr(10)=+cr(3)*cr(1).^2;
v(11)=v(2);t0(11)=t0(3)+3*t0(1);cr(11)=-cr(3)*cr(1).^3;
%v(12)=v(1);t0(12)=1.8;cr(12)=+cr(1)^9;
%v(13)=v(1);t0(13)=2.0;cr(13)=+cr(1)^10;
end
nr=max(size(v));

if optionw=='rick'
w=myricker(nw,f,dt);
figure,plot(w)
elseif optionw=='minr'
   load minricke.dat;
w=minricke;   
end


nt=round(tf/dt)+1;
xx=zeros(2*nt,2*nh);

for rr=1:nr;
for hh=1:nh
   h=dh*(hh-1)+h_near;
   t=sqrt(t0(rr).^2+(h/v(rr)).^2);
   tt=round(t./dt)+1;
   %if hh==23 tt, end
   xx(tt,hh)=xx(tt,hh)+cr(rr);
end;
end;
t=0:dt:tf;
tt=t./dt+1;
tt=round(tt);
for rr=1:nr;
   for hh=1:nh
      xxtemp=conv(xx(tt,hh),w);
      xxc(1:nt,hh)=xxtemp(fix(nw/2)+1:nt+fix(nw/2));
	end;
end;

if optionn=='y';
   noise=rand(size(xxc));
   xxc=xxc+noise*0.000;
end;

xxtemp=xx./max(max(abs(xx)));
h=(0:nh-1)*dh;
figure,wigb(real(xx(1:nt,1:nh)),1,h,t);
%figure;plotrace(xxtemp,dh,nh,h_near)
figure,wigb(real(xxc(1:nt,1:nh)),1,h,t);

xxtemp=xxc./max(max(abs(xxc)));
%figure;plotrace(xxtemp(1:nt,:),dh,nh,h_near)
ylabel('time (sec)');xlabel('offset (m)');
title('Synthetic Seismic CDP');
if options=='y' save c:\daniel\danilo.mat xxc /ascii;end;
