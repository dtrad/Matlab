clear;
V(1)=1000;
V(2)=1500;
V(3)=2000;
V(4)=3000;

% densities
r0=0;
r(1)=1000;
r(2)=1300;
r(3)=2000;
r(4)=3000;

% two way travel distances;
z(1)=200;
z(2)=300;
z(3)=250;
%h(4)=600;

NT=512; % Number of time
NP=80;  % Number of p traces.
NH=40;
dh=10;h=(0:(NH-1)).*dh;
dp=1/(NP*V(1)); % Delta p
dt=0.004;

nl=length(V);
u=1./V;
m=1:NP;
p=(m-1)*dp;
j=1:nl;
um=u(:)*ones(1,length(p));
pm=ones(length(u),1)*p(:).';
q(j,m)=(um.^2-pm.^2).^0.5;
R=zeros(NT,NH);
for n=1:nl;
	for m=1:NP;
   	j=1:n-1;
   	tau(n,m)=2*(z(j)*q(j,m));
      time=t2index(tau(n,m),dt);
      R(time,m)=1;
   end
end
R(1,:)=0;
wav=rickerm(50,dt);
seis_shape(R);
figure,
wigb(R);
for m=1:NP;Rw(:,m)=convlim(R(:,m),wav,NT);end

xt=taup2xt(Rw,'linear',h,p,dt);
figure,
wigb(xt),
% First response R(nl)=0
