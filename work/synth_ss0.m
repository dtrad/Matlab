% model
clear;
close all
options='n'; %save
options2='y'; %save
optionm='y'; %multiples
%optionw='minr'; % wavelet
optionw='rick';
optionn='y'; % noise;

tf=2.2;
dt=0.004;
dh=25;
nh=50;
h_near=-1250;
nw=50;
f=100;
aten=-1.2;
aten=0;

v(1)=2500;t0(1)=0.3;cr(1)=0.5;
v(2)=3300;t0(2)=1.7;cr(2)=(1-cr(1)^2)*0.5;%cr(2)=0;
v(3)=v(2);t0(3)=0.8;cr(3)=0.;
v(4)=v(2);t0(4)=1.27;cr(4)=0.0;

% multiples
if optionm=='y'
v(5)=v(1);t0(5)=2*t0(1);cr(5)=-cr(1)^2;
v(6)=v(1);t0(6)=3*t0(1);cr(6)=+cr(1)^3;
v(7)=v(1);t0(7)=4*t0(1);cr(7)=-cr(1)^4;
v(8)=v(1);t0(8)=5*t0(1);cr(8)=+cr(1)^5;
v(9)=v(2);t0(9)=2*t0(2)+t0(1);cr(9)=-cr(2).^2*cr(1);
v(10)=v(2);t0(10)=t0(2)+2*t0(1);cr(10)=-cr(2)*cr(1).^2;
%v(11)=v(2);t0(11)=t0(3)+3*t0(1);cr(11)=-cr(3)*cr(1).^3;
%v(12)=v(1);t0(12)=1.8;cr(12)=+cr(1)^9;
%v(13)=v(1);t0(13)=2.0;cr(13)=+cr(1)^10;
end
nr=max(size(v));

if optionw=='rick'
w=rickerm(f,dt);%w=padzeros(w,nw);
figure,plot(w)
elseif optionw=='minr'
   load minricke.dat;
w=minricke;   
end


nt=round(tf/dt)+1;
xx=zeros(2*nt,2*nh);

for rr=1:nr;
   for hh=1:nh
      h(hh)=dh*(hh-1)+h_near;
      %avo=abs(cos((h-h_near)./((nh-1)*dh)*pi));
      %t=sqrt(t0(rr).^2+(h/v(rr)).^2);
      t=t0(rr)+(abs(h(hh))/v(rr));
      tt=(t./dt)+1;
      tt1=tt-floor(tt);
      tt2=1-tt1;
   	%if hh==23 tt, end
      %xx(tt,hh)=(xx(tt,hh)+cr(rr)).*avo;
     	xx(floor(tt),hh)=(xx(floor(tt),hh)+tt2*cr(rr));
     	xx(ceil(tt),hh)=(xx(ceil(tt),hh)+tt1*cr(rr));
	end;
end;
tt=1:nt;
for rr=1:nr;
   for hh=1:nh
      xxtemp=conv(xx(tt,hh),w);
      xxc(1:nt,hh)=xxtemp(1:nt);
      %xxc(1:nt,hh)=xxtemp(fix(nw/2)+1:nt+fix(nw/2));
	end;
end;

if optionn=='y';
   noise=randn(size(xxc));
   xxc=xxc+(noise-0.5)*0.001;
end;

xxtemp=xx./max(max(abs(xx)));
t=(0:nt-1)*dt;ntlim=ceil(nt/2);
figure,wigb(real(xx(1:ntlim,1:nh)),1,h,t(1:ntlim));
%figure;plotrace(xxtemp,dh,nh,h_near)
figure,wigb(real(xxc(1:ntlim,1:nh)),1,h,t(1:ntlim));

xxtemp=xxc./max(max(abs(xxc)));
%figure;plotrace(xxtemp(1:nt,:),dh,nh,h_near)
ylabel('time (sec)');xlabel('offset (m)');
title('Synthetic Seismic CDP');
if options=='y' save c:\daniel\thesis\xxc_SS.mat xxc xx dh dt h_near;end;
if options2=='y'; 
   xxc=xxc(1:512,1:50);
	xxc=xxc(:);
	h=h(:);
   save c:\daniel\radon\ss0.txt xxc /ascii;
   save c:\daniel\radon\ss0.off h /ascii;
end;
