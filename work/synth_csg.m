% model
clear;
close all;
optionm='y'; %multiples
%optionw='minr'; % wavelet
optionw='rick';
optionn='n'; % noise;
options='y';

tf=2.2;
dt=0.004;
dh=10;
ngeop=64;
nshot=64;
xs=1:nshot;xs=(xs-1)*dh;
xg=1:ngeop;xg=(xg-1)*dh;
nw=50;
f=70;
aten=-1.2;
aten=0;

v(1)=1000;t0(1)=0.4;cr(1)=0.3;

% multiples
if optionm=='y'
v(2)=v(1);t0(2)=2*t0(1);cr(2)=-cr(1)^2;
v(3)=v(1);t0(3)=3*t0(1);cr(3)=+cr(1)^3;
v(4)=v(1);t0(4)=4*t0(1);cr(4)=-cr(1)^4;
v(5)=v(1);t0(5)=5*t0(1);cr(5)=+cr(1)^5;
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
hh=1:ngeop;

xx=zeros(nt,ngeop,nshot);
for ns=1:nshot;
   for rr=1:nr;
      for hh=1:ngeop;
      	h=xg(hh)-xs(ns);
   		t=sqrt(t0(rr)^2+(h/v(rr)).^2);
   		tt=round(t./dt)+1;
      	xx(tt,hh,ns)=(xx(tt,hh,ns)+cr(rr));
      end
	end;
end   

[nt,nr,nshot]=size(xx);

for ns=1:nshot
   for rr=1:nr
      xx(:,rr,ns)=convlim(xx(:,rr,ns),w,nt);
   end   
end   
   


t=0:dt:tf;
tt=t./dt+1;
tt=round(tt);

%for hh=1:nh;xxtemp=conv(xx(tt,hh),w);xxc(1:nt,hh)=xxtemp(1:nt);end;

if optionn=='y';
   noise=rand(size(xxc));
   xx=xx+noise*0.000;
end;


for ns=1:3:ns
   h=xg-xs(ns);
	figure,wigb(real(xx(1:nt,1:ngeop,ns)),1,h,t);
end

ylabel('time (sec)');xlabel('offset (m)');
title('Synthetic Seismic CDP');

if options=='y' save c:\daniel\thesis\shot_gathers.mat xx;end;
