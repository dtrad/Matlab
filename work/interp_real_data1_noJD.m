% Interpolation for real data
% First time read the data from segy file
pass=2;
if pass=='1'
[D,H]=readsegy('some_cdps.segy');
[D,H]=ssort(D,H,'offset');
D=seis_shape(D);
Ro=D(1:512,50:70);
h0=[H.offset];
h0=h0(20:100);
save c:\daniel\thesis\real_data.mat Ro h0;
else
   load c:\daniel\thesis\real_data;
end

z=seis_shape(Ro(1:512,:));
%z=mute_matrix(z);
z=normalize(z);
[nt,nh]=size(z);
dh=20;

%[h1,gap]=h1_without_gaps(h0,dh);

h1=min(h0):dh:max(h0);

dt=0.004;
pnorm=1;
option='Yilmaz ';
vel_min=800;
vel_max=5000;
np=80;dp=1/vel_min/np;
niter=15;
hyperparameter=1e-4;
optionmute='y';
freqint=[2 255];

[u,haxis,ttaxis,vr,ppaxis]=interp_taup1_noJD(z,h0,h1,dt,pnorm,option,vel_min,np,niter,hyperparameter,vel_max,optionmute,freqint);
[vel,pr]=alfa2vel(ppaxis,[0 dp]);
pr=padzeros(pr,length(ppaxis));
figure,wigb(vr,1,ppaxis,ttaxis);
[vel,pr]=alfa2vel(ppaxis,[0 dp]);
pr=padzeros(pr,length(ppaxis));
subplot(211),semilogy(vel);
vel(1:10)=vel(1:10)*1e-5;
subplot(211),semilogy(vel);title('Vel');xlabel('#trace')
subplot(212),semilogy(ppaxis);title('alfa');title('alfa');xlabel('#trace')




 