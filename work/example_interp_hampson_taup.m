load output;
z=seis_shape(xxc(1:512,1:30));
[nt,nh]=size(z);

temp1=(0:9)*50;
temp2=(10:19)*50; % Gap
temp3=(20:29)*50;

h0=[temp1 temp3];
h1=[temp1 temp2 temp3];

select=[ones(1,10) zeros(1,10) ones(1,10)];
x=shrinkt2(z,select);

[u,haxis,ttaxis,vr]=interp_taup2(x,h0,h1,[],1,'Hampson',1000,40,10,1e-3);
