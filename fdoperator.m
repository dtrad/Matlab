function fdoperator(v)
if (nargin < 1) v=0.1;end
nt=4;
nz=5;
loc=round(nz/2+0.5);
loc=nz-1;
d=zeros(nt,1);
m=zeros(nz,1);
m(loc)=1;
mm=m*ones(1,nt);
s0=zeros(nt,nz);
s1=s0;
s2=s0;
s3=s0;
s4=s0;
s0(1,1)=1;
s1(2,1)=1;
s2(3,1)=1;
s3(4,1)=1;

[s0;s1;s2;s3]

a=v*2*ones(nz,1);

t1=diag(a,-1);t1=t1(1:nz,1:nz);
t2=diag(2*(1-a));t2=t2(1:nz,1:nz);
t3=diag(a,1);t3=t3(1:nz,1:nz);
T=t1+t2+t3;
%T=diag(2*(1-a));

I=eye(nz);
O=zeros(nz,nz);
FD3=[I O O O;O I O O; O O I O;O -I T O];
FD2=[I O O O;O I O O;-I T O O;O  O O O];
FD1=[I O O O;T O O O; O O O O;O  O O O];
FD0=[I;O;O;O];

S=[s0 s1 s2 s3];
OP1=FD3*FD2*FD1*FD0*mm;
OP=S*OP1
figure(1);
subplot(221);imagesc(S);colorbar;
subplot(222);imagesc(OP1);colorbar;
subplot(223);plot(OP);title('data');
subplot(224);plot(m,'o');title('reflectivity')
figure(gcf)
return;



