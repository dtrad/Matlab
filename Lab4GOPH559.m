%GOPH 559 lab4



dt=0.001;
v1=1500;h1=100;r1=1;
v2=2600;h2=250;r2=2.5;
v3=2000;h3=200;r3=2;
v4=3000;h4=200;r4=3;
re1=refl(v1,v2,r1,r2,h1,dt);
re2=refl(v2,v3,r2,r3,h2,dt);
re3=refl(v3,v4,r3,r4,h3,dt);
rtotal=[re1;re2;re3;zeros(floor(h4/v4/dt),1)]

w=ricker(dt);
t=conv(rtotal,w,'same');
subplot(311);plot(w);
subplot(312);plot(rtotal);
subplot(313);plot(t);

figure(gcf);

function [r]=refl(v1,v2,r1,r2,h1,dt)
c=(v2*r2-v1*r1)/(v2*r2+v1*r1);
nzeros=floor((h1/v1)/dt);
r=zeros(nzeros,1);
r(end)=c;
end