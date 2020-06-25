function xnoah=noahdecon(R,wav)
% xnoah=noahdecon(R)
% Noah deconvolution:
% Given the trace R xnoah is th deconvolution of 1+R from R
% xnoah=R/(1+R)
% Since xnoah is the trace which would had been recorded is the free 
% surface would not exist, xnoah will be free of surface multiples
% and only peg legs remain.
% Reference; Claerbout, Geophysics. vol 41 (august 1996)
% Daniel Trad- UBC- 23/07/98
option='n';
if nargin==2 option='w';elseif nargin==1 option='s';end
if option=='s'
lx=length(R);
R(1)=0;
Rplus1=-R;Rplus1(1)=1;x=zeros(lx,1);x(1)=1;
xnoah=filter(1,Rplus1,R);
xnoah=xnoah(1:lx);
elseif option=='w'
   lx=length(R);
   %ref=zeros(lx,1);ref(1)=1;ref=convlim(ref,wav,lx);R=R-ref;
   R_decon=deconv(R,wav);
   Rplus1=-R_decon;Rplus1(1)=1;
   x=zeros(lx,1);x(1)=1;
	xnoah=filter(1,Rplus1,R);
	xnoah=xnoah(1:lx);
elseif option=='n'
   display('wrong number of arguments')
end   
