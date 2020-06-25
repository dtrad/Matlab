% Tomography problem: Geop507- Exercise 7
% Solving and Overdetermined System of Equations
% All the cells of the model are arranged in one single vector
% sorted by rows (i.e., one row after another)
% x is the horizontal dimension and y is the vertical
% kx and ky are the sizes, dx dy the intervals
% 
% This version is a modification of geop507c. The same
% inversion is carried out but setting the problem as overdetermined.
% This is, Ax=b means 
%		/						\
% 		|WdG					|
% 	A=	|sqrt(alfas) Ws	|
%		|sqrt(alfax) Wx	|
%		|sqrt(alfaz) Wz	|
%		\						/
% Daniel Trad


close all
clear
% Geometry 0
yl=5:10:95;
xl=0;
yr=5:10:95;
xr=100;
vel=1000;
vel0=1000;  % Reference model
vel2=2000;
alfas=0.01;
alfax=1;
alfay=1;
sref=1;
OPTIONSVD='n';
% Noise 
randn('seed',1);
stdnoise=2e-3;

% Define model
xmin=50;
xmax=80;
zmin=50;
zmax=80;

%Geometry 0
%yl=[10,30,50,70,90];
%xl=0;
%yr=[10,30,50,70,90];
%xr=100;

dx=10;
dy=10;
kx=10;
ky=10;

% Construction of A using the raypath
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

slowness=1/vel;
s=slowness*ones(kx*ky,1);
% Reference model
slowness0=1/vel0;
m0=sref*slowness0*ones(kx*ky,1);

%s=reshape(s,kx,ky);
s=mat(s,kx);
s((zmin/dy+1):(zmax/dy),(xmin/dx+1):(xmax/dx))=1/vel2;

figure,
imagesc(s);colorbar;title('model');      

%s=s(:);
s=vc(s);
to=A*s;
noise=stdnoise*randn(size(to));
t=to+noise;
tindex=1:length(t);
figure,plot(tindex,to,tindex,t,'o'),title('exact data ---, noise data'); 


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
   x(round(t(ii)/dt+1),ii)=1;
end

figure,wigb(x,1,xaxis,taxis);

rr=rank(A)
display(sprintf('Rank of A=%d\n',rr));
% Data Covariance Matrix
Wd=diag(1./stdnoise);

if OPTIONSVD=='y'
   % Adjoint
   sa=A'*t;
	scale=median(diag(A'*A));
	sa=sa./scale;

	figure,
	subplot(211),plot(s),title('model')
	subplot(212),plot(sa),title('adjoint')
   
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
	[(t-tp)'*(t-tp),(t-tpc)'*(t-tpc),(t-tpa)'*(t-tpa)]

	display('Norms: True model//Constr Model//Adjoint model'); 
	[s(:).'*s(:),mc(:).'*mc(:),sa(:).'*sa(:)]

	mpinv=pinv(A)*t;

	s=mat(s,kx);
	sa=mat(sa,kx);
	mc=mat(mc,kx);
	mpinv=mat(mpinv,kx);

	figure,
	subplot(221); imagesc(s);colorbar;title('model');      
	subplot(222); imagesc(mc);colorbar,title('constructed model');      
	subplot(223); imagesc(sa);colorbar,title('adjoint model');      
	subplot(224); imagesc(mpinv);colorbar,title('Pseudo inverse model');      

	figure,
	subplot(211);plot(U(:,rr+1:length(s(:))));title('Eigenvector for the inactivated data space')
	subplot(212);plot(V(:,rr+1:length(s(:))));title('Eigenvector for the inactivated model space')
   
	ttaux=Wd*(t-tpc);
	phi=ttaux(:).'*ttaux;	
	display(sprintf('True misfit=%e\n',phi)); 
end


% dx and dy are constants, otherwise we need for loops
dxdy=sqrt(dx*dy);
Wms=dxdy*eye(length(s(:)));

% Build Wmx= dm/dx
Wmx=eye(size(Wms));
[nr nc]=size(Wmx);
for ii=1:nr-1
   if (mod(ii,kx)~=0) Wmx(ii,ii+1)=-1;
   else Wmx(ii,ii)=0;
	end  
end
Wmx=dxdy*Wmx;

% Build Wmy= dm/dy
Wmy=eye(size(Wms));
[nr nc]=size(Wmy);
for ii=1:nr-kx
   Wmy(ii,ii+kx)=-1;
end
for ii=nr-kx+1:nr
   Wmy(ii,ii)=0;
end
Wmy=dxdy*Wmy;

% Total matrix for the left hand side (F stands for forward problem)
F=[Wd*A;sqrt(alfas)*Wms;sqrt(alfax)*Wmx;sqrt(alfay)*Wmy];
% Right hand side
b=[t;sqrt(alfas)*Wms*m0;sqrt(alfax)*Wmx*m0;sqrt(alfay)*Wmy*m0];
%Wm=alfas*Wms+alfax*Wmx+alfay*Wmy;
%Wm=eye(length(mc(:)));
%norm(Wm*vc(mc));
beta0=1;beta1=8;
beta=logspace(beta0,beta1,9);
%beta=[1e-3 1e-1 10 1e2 1e3 1e4 1e5 1e6 1e7];
figure,
for ii=1:9
   mc=inv(F'*F+beta(ii)*eye(length(s)))*F'*Wd*b;
   %mc=minphi(A,Wm,Wd,t,beta(ii),m0);
   tp=A*mc;
   phim(ii)=alfas*(Wms*(mc-m0)).'*(Wms*(mc-m0))+...
      alfax*(Wmx*(mc-m0)).'*(Wmx*(mc-m0))+...
      alfay*(Wmy*(mc-m0)).'*(Wmy*(mc-m0));
   phid(ii)=(Wd*(t-tp)).'*(Wd*(t-tp));
   phit(ii)=phid(ii)+beta(ii)*phim(ii);
   subplot(330+ii),imagesc(mat(mc,kx));colorbar;
end
figure
subplot(221);loglog(phim,phid,'o'),title('L curve');
subplot(222);loglog(beta,phid,'o'),title('\phi_d');
subplot(223);loglog(beta,phim,'o'),title('\phi_m');
subplot(224);loglog(beta,phit,'o'),title('\phi');

% Estimation of Beta
chi=length(s(:));
chiup=chi*1.1;
chidown=chi*0.9;

xx=logspace(beta0,beta1,100);
yy=spline(beta,phid,xx);

iless=find(yy<=chiup);
ibeta=find(yy(iless)>=chidown);
betaf=mean(xx(ibeta))

mc=inv(F'*F+betaf*eye(length(s)))*F'*Wd*b;
tp=A*mc;
phimf=alfas*(Wms*(mc-m0)).'*(Wms*(mc-m0))+...
      alfax*(Wmx*(mc-m0)).'*(Wmx*(mc-m0))+...
      alfay*(Wmy*(mc-m0)).'*(Wmy*(mc-m0));
phidf=(Wd*(t-tp)).'*(Wd*(t-tp));
phif=phidf+betaf*phimf;
display(sprintf('Final misfit=%e\n',phidf)); 
display(sprintf('final norm=%e\n',phimf));
display(sprintf('cost function=%e\n',phif));
figure,
subplot(211);plot(t,'o'),title('Data')
subplot(212);plot(tp,'o'),title('Predicted by the constructed model')
figure,imagesc(mat(mc,kx));colorbar;

figure,
subplot(211),plot(svd(F))
subplot(212),plot(svd(A*Wd'*Wd*A+alfas*Wms'*Wms+alfax*Wmx'*Wmx+alfay*Wmy'*Wmy))

