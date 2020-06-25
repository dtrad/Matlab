% Tomography problem: Geop507
% Daniel Trad
	
close all
clear
% Geometry 0
yl=[10,30,50,70,90];
xl=0;
yr=[10,30,50,70,90];
xr=100;

% Geometry 0
yl=[10,30,50,70,90];
xl=0;
yr=[10,30,50,70,90];
xr=100;


dx=20;
dy=20;
kx=5;
ky=5;

for ii=1:length(yl)
	for jj=1:length(yr)
	%hold on
	[v]=ray(xl,yl(ii),xr,yr(jj),dx,dy,kx,ky);
	%pause
	if (ii==1)&(jj==1) A=v;
	else A=[A;v];
	end

	end
end
%hold off

figure,
imagesc(A)      
colorbar('ver')
colormap(jet)

vel=1000;
slowness=1/vel;
s=slowness*ones(kx*ky,1);
%s=reshape(s,kx,ky);
s=mat(s,5);

% Define model
xmin=40;
xmax=80;
zmin=40;
zmax=100;

s((zmin/dy+1):(zmax/dy),(xmin/dx+1):(xmax/dx))=1/2000;

figure,
imagesc(s);colorbar
title('model');      

%s=s(:);
s=vc(s);
t=A*s;
% Seismogram    (wiggle traces)
tmax=max(t)+.1;
dt=tmax/100;;
nt=round(tmax/dt);
nh=length(t);

taxis=0:nt-1;
taxis=taxis*dt;
xaxis=1:nh;

x=zeros(nt,length(t));
for ii=1:nh
   x(round(t(ii)/dt),ii)=1;
end

figure,wigb(x,1,xaxis,taxis);

% Adjoint

sa=A'*t;

scale=median(diag(A'*A));

sa=sa./scale;

figure,
subplot(211),plot(s),title('model')
subplot(212),plot(sa),title('adjoint')


display('Rank of A')
rr=rank(A)

[U,S,V]=svd(A,0);
figure,
plot(diag(S));title('Singular Values of A')


mc=V(:,1:rr)*((U(:,1:rr)'*t)./diag(S(1:rr,1:rr)));
size(mc)
mc=[mc;zeros(length(s)-length(mc),1)];

figure,
subplot(311);plot(s,'o'),title('model')
subplot(312);plot(mc,'o'),title('constructed model')
subplot(313),plot(sa,'o'),title('adjoint')

tp=A*s;
tpc=A*mc;
tpa=A*sa;

figure,
subplot(311);plot(tp,'o'),title('Predicted by the model')
subplot(312);plot(tpc,'o'),title('Predicted by the constructed model')
subplot(313),plot(tpa,'o'),title('Predicted by the adjoint')

format short e
display('Misfits True model//Constr Model//Adjoint model'); 
[norm(t-tp),norm(t-tpc),norm(t-tpa)]

display('Norms: True model//Constr Model//Adjoint model'); 
[norm(s(:)),norm(mc(:)),norm(sa(:))]


mpinv=pinv(A)*t;

s=mat(s,5);
sa=mat(sa,5);
mc=mat(mc,5);
mpinv=mat(mpinv,5);

figure,
subplot(221); imagesc(s);colorbar;title('model');      
subplot(222); imagesc(mc);colorbar,title('constructed model');      
subplot(223); imagesc(sa);colorbar,title('adjoint model');      
subplot(224); imagesc(mpinv);colorbar,title('Pseudo inverse model');      

figure,
subplot(211);plot(U(:,rr+1:length(s(:))));title('Eigenvector for the inactivated data space')
subplot(212);plot(V(:,rr+1:length(s(:))));title('Eigenvector for the inactivated model space')








