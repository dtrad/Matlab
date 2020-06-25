clear all;
close all;
%velocity model
zp=0:10:2500;
vp(1:100)=1860+.6*zp((1:100)); 
vp(101:155)=3100+0.4*(zp(101:155)-1000);
vp(156:251)=2750+0.7*(zp(156:251)-1550);
vs=.5*vp;zs=zp;
plot(vp,zp);flipy; 
xlabel('m/s'); ylabel('meters'); title('v(z) model');

%setting parameters
zsrc=0;zrec=0;zd=2500;%source receiver and reflector depths
xoff=10:100:3010;caprad=10;itermax=4;%offsets,cap radius,and max iter
pfan=-1;optflag=1;pflag=1;dflag=2;%default ray fan,and variousflags

%create P-P reflection
figure;subplot(2,1,1);flipy;
[t,p1]=traceray_pp(vp,zp,zsrc,zrec,zd,xoff,caprad,pfan,itermax,optflag,pflag,dflag);
title(['P-P Reflection'])
grid;xlabel('meters');ylabel('meters');

%calculate and plot the angle of incidence on the reflector versus offset
angle = asin(1860*p1)*(180/pi);
subplot(2,1,2);plot(xoff,angle);grid;
xlabel('meters');ylabel('degrees')

%P-S reflection
figure;subplot(2,1,1);flipy;
[t,p2]=traceray_ps(vp,zp,vs,zs,zsrc,zrec,zd,xoff,caprad,pfan,itermax,optflag,pflag,dflag);
title(['P-S Reflection'])
grid;xlabel('meters');ylabel('meters ');

%calculate and plot the angle of incidence on the reflector versus offset
angle = asin(1860*p2)*(180/pi);
subplot(2,1,2);plot(xoff,angle);grid;
xlabel('meters');ylabel('degrees')