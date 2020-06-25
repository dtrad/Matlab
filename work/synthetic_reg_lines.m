% model
clear;
for ii=1:15,close,end;
options='y'; %save
optionm='n'; %multiples
%optionw='minr'; % wavelet
optionw='rick';
optionn='y'; % noise;

tf=2.2;
dt=0.004;
dh=25;
nh=70;
h_near=0;
nw=50;
f=25;
aten=-1.2;

v(1)=2000;t0(1)=0.4;cr(1)=0.3;
v(2)=-1300;t0(2)=0.8;cr(2)=0.3;
%v(3)=1500;t0(3)=1;cr(3)=-0.2;
%v(4)=-1000;t0(4)=1.5;cr(4)=0.2;

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

%generate irregular offset.
for hh=1:nh
   if (hh==1) h(1)=h_near;
   else h(hh)=h(hh-1)+(1+0.*rand(1))*dh;
   end
end
   
%hh=1:nh;h=dh*(hh-1)+h_near;
   
for rr=1:nr;
for hh=1:nh

   %t=sqrt(t0(rr).^2+(h(hh)/v(rr)).^2);
   t=t0(rr)+h(hh)/v(rr);
   tt=round(t./dt)+1;
   %if hh==23 tt, end
   if (t>=0) 
       xx(tt,hh)=xx(tt,hh)+cr(rr);
   end
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
%h=(0:nh-1)*dh;
figure,wigb(real(xx(1:nt,1:nh)),1,h,t);
%figure;plotrace(xxtemp,dh,nh,h_near)
figure,wigb(real(xxc(1:nt,1:nh)),1,h,t);

xxtemp=xxc./max(max(abs(xxc)));

data = xxtemp;clear xxtemp;
%figure;plotrace(xxtemp(1:nt,:),dh,nh,h_near)
ylabel('time (sec)');xlabel('offset (m)');
title('Synthetic Seismic CDP');
if options=='y' save c:\dtrad\MATLAB7\work\lineregnoalias.mat data h t;end;
