load c:\albert\zz;
z=seis_shape(z);
[nt,nh]=size(z);

temp1=(0:9)*25;
temp2=(10:19)*25; % Gap
temp3=(20:31)*25;

h0=[temp1 temp3];
h1=[temp1 temp2 temp3];

select=[ones(1,10) zeros(1,10) ones(1,12)];
x=shrinkt2(z,select);

[u,haxis,ttaxis,vr]=interp_taup2(x,h0,h1,[],1,'linear ',1000,40,10,1e-3);
