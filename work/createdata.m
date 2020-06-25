function [data,datart,h,q,t]=createdata(nt,nh,nq,fp,dt,dh,dq,iq,it)
% Function to create parabolic or linear data from a Radon space.
%
%            data=createdata(nt,nq,fp,dt,iq,it)
%
% This is usuful for testing RT algorithms. 
% input 
%      nt
%      nq
%      fp 
%      dt 
%      iq 
%      it
% Output 
%       data (t,x)
%       datart (tau,q)

if (nargin<1) nt=512;end
if (nargin<2) nh=60;end
if (nargin<3) nq=40;end
if (nargin<4) fp=25;end
if (nargin<5) dt=0.004;end
if (nargin<5) dh=10;end
if (nargin<5) dq=5e-7;end
if (nargin<6) iq=[round(nq/2)-1,round(nq/2)+1];end
if (nargin<7) it=ones(2,1)*round(nt/3);end

ne=length(iq);
m=zeros(nt,nq);

for i=1:ne
  m(it(i),iq(i))=1;
end

w=rickerm(fp,dt);
nw=length(w);
data=zeros(nt+nw-1,nq);

for i=1:nq
  mm=m(:,i);
  datart(:,i)=conv(mm(:),w(:));
end

datart=datart(1:nt,1:nq);

q=0:nq-1;
q=q*dq;
h=0:nh-1;
h=h*dh;h=h-h(round(nh/2));
t=0:nt-1;
t=t*dt;

data=radonL(datart,q,dt,h);

