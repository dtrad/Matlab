% model
clear;
for ii=1:15,close,end;
options='y';

tf=2;
dt=0.004;
dh=50;
nh=70;
h_near=0;
nw=51;
f=40;
aten=-1.2;


v(1)=3000;t0(1)=0.2;cr(1)=0.5;
v(2)=3300;t0(2)=0.43;cr(2)=0.0;
v(3)=v(2);t0(3)=0.93;cr(3)=0.4;
v(4)=v(2);t0(4)=1.27;cr(4)=0.0;

% multiples
v(5)=v(1);t0(5)=0.4;cr(5)=-cr(1)^2;
v(6)=v(1);t0(6)=0.6;cr(6)=+cr(1)^3;
v(7)=v(1);t0(7)=0.8;cr(7)=-cr(1)^4;
v(8)=v(1);t0(8)=1.0;cr(8)=+cr(1)^5;
v(9)=v(1);t0(9)=1.2;cr(9)=-cr(1)^6;
v(10)=v(1);t0(10)=1.4;cr(10)=+cr(1)^7;
v(11)=v(1);t0(11)=1.6;cr(11)=-cr(1)^8;
v(12)=v(1);t0(12)=1.8;cr(12)=+cr(1)^9;
%v(13)=v(1);t0(13)=2.0;cr(13)=+cr(1)^10;

nr=max(size(v));


w=ricker(nw,f,dt);
figure,plot(w)

nt=round(tf/dt)+1;
xx=zeros(nt,nh);

for rr=1:nr;
for hh=1:nh
   h=dh*(hh-1)+h_near;
   t=sqrt(t0(rr).^2+(h/v(rr)).^2);
   tt=round(t./dt)+1;
   xx(tt,hh)=cr(rr);
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
xxtemp=xx./max(max(abs(xx)));
figure;plotrace(xxtemp,dh,nh,h_near)
xxtemp=xxc./max(max(abs(xxc)));
figure;plotrace(xxtemp(1:nt,:),dh,nh,h_near)
xlabel('time (sec)');ylabel('offset (m)');
if options=='y' save c:\daniel\taup521\output.mat xxc ;end;
