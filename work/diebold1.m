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
z(1)=400;
z(2)=600;
z(3)=450;
%h(4)=600;

NT=1024; % Number of time
NP=80;  % Number of p traces.
NH=40;
dh=10;
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
for xx=1:NH
x=(xx-1)*dh;
for n=1:nl;
	for m=1:NP;
   	j=1:n-1;
   	tau(n,m)=2*(z(j)*q(j,m));
   	t(n,m)=p(m)*x+tau(n,m);
   end
   tx(n,xx)=sum(t(n,:));
   time=t2index(tx(n,xx),dt);
	R(time,xx)=1;
end
end
seis_shape(R);
wigb(R);

% First response R(nl)=0
%wav=rickerm(50,0.004);
