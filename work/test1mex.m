mex -g c:\daniel\cpp\forwardmex4.cpp
UH=[(2+i),(3-i); (3-i),4; (1+i) 2; 2 3;  3 5; (2-i),3];
VLS=zeros(6,2);VLS(1,1)=1+i;
tol=0.001;
NITER=2;
tolmodel=0.001;
pnorm=1;
sigmap=1;
Cn=ones(6,2);
w=[1 1 2 3 4 5];
h0=[1 1];
nt=2;
Power=1;
WV=ones(2,2);
WU0=ones(2,2);
W=ones(2,2);
np=2;
alfa=[1 1];
freqint=[0 5];

[V]=forwardmex4(UH,VLS,tol,NITER,tolmodel,pnorm,sigmap,Cn,w,h0,nt,Power,WV,WU0,W,np,alfa,freqint)