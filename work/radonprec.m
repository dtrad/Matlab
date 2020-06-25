function [V,DR]=radonprec(D,w,h,p,method,vtrue,rmethod)
% Test different algorithms for precomputing the least square
% operator for the Radon transforms with irregular p sampling
%
% Given the data D computes the Radon transform
% for every frequency w(i)
% h is the offset axis and p the radon axis.
% vtrue is the model from which the data were generated.
% If available the true model is used to compute the model covariance
% Qp, allowing a high resolution RT
% method is one of the following
% 1  Conjugate gradients Square with Tikhonov
% 2  Conjugate gradients Least Squares (niter gives reg)
% 3  Weighted Conjugate gradients Least Squares   
% 4  wlsqr (Weigthed BD) with GCV.
% 5  Truncated SVD
% 6  Conjugate gradients Least Squares
% 7  Truncated SVD and GCV.
% 8  Tikhonov with SVD with GCV.
% 9  Truncated GSVD with GCV.
% 10 Tikhonov with GSVD ans GCV.
% 11  Weighted Conjugate gradients Least Squares with eigenvalue plot
% 12 Weighted Conjugate gradients Least Squares with regularization demo
% 
% Daniel Trad-- September - 02
if nargin<7|isempty(rmethod) rmethod='LRT';end  

global stpc  reorth eps2 step;
stpc
reorth
eps2
step
rmethod

nh=length(h);
dh(1)=h(2);ii=2:nh;
dh(ii)=h(ii)-h(ii-1);
np=length(p);
V=zeros(1,np);
DR=zeros(1,nh);
nf=length(w);

% Preallocate memory
FH3=zeros(np,nh,nf);
F3=zeros(nh,np,nf);

Cm=zeros(nf,np);
Qp3=zeros(nf,np);
R3=zeros(np,np,nf);
LS3=zeros(np,nh,nf);

%Qp=ones(1,np);
ke=zeros(1,length(w));
index=0;
Vtrue=fft(vtrue);

% if vtrue is not available use eps2=1;
if (max(max(abs(vtrue)))<1e-10) eps2=1;end
condn=zeros(nf,1); 
II=eye(np);

% Estimate the variance of the noise using a high frequency with low 
% signal noise ratio
d=D(nf-10,:).';
eps1=real(d'*d)/nh;  % This is the variance of the noise
message=sprintf('Variance of the noise eps1=%e\n',eps1);
display(message);

% Precalculation of LS operators
sigma=1/sqrt(eps1);
Cd=sqrt(eps1)*ones(size(d));
tic
for f=2:nf-5,
  % Construct and store the general operators required by all methods 
  [FH,F]=radonmat2(w(f),h,p,rmethod);
   FH3(:,:,f)=FH;
   F3(:,:,f)=F;
   %condn(f)=cond(FH*F+1e-11*II);
   %if (f==30) keyboard;end
   %FH=F';
   %FH=F'; %'
   Cm(f,:)=(eps2+abs(Vtrue(f,:)));
   Qp3(f,:)=sqrt(1./(eps2+abs(Vtrue(f,:)))+1e-7);



   if (method==1)
     % precalculation with Cholesky operator
     % force the diagonal to be real for cholesky
     FHF=FH*F;
     for ip=1:np,
       FHF(ip,ip)=real(FHF(ip,ip));
     end
     

     [R]=chol(FHF+diag(1./Cm(f,:)));
     R3(:,:,f)=R;
   elseif (method==2)
     % precalculate the whole LS operator    
     FHF=FH*F;
     [LS]=(inv(FHF+diag(1./Cm(f,:))))*FH;
     LS3(:,:,f)=LS;
   end  
end   
toc
display('end of precomputation');

% Loop using the precalculated LS operators
for icmp=1:1, % The same calculation is done for several CMPS
  display('beginning of cmp calculation');
  tic
  V=zeros(1,np);
  for f=2:nf-30,   
    %if f==nf-27 [R]=resolution(FH,F,Cd,Cm,rmethod,p);end
    d=D(f,:).';  %'
    FH=FH3(:,:,f);
    F=F3(:,:,f);
    Qp=Qp3(f,:);
    if (method==0)
      [v]=FH*d;
    elseif (method==1)
      %[v,flag]=cgs(FH*F+0.01*eye(np),FH*d);
      R=R3(:,:,f);
      va=FH*d;
      v= R\(R'\va);
    elseif (method==2)
      %niter=round(f/nf*np)
      
      niter=10;
      [v]=LS3(:,:,f)*d;
   elseif (method==3)
      WW=diag(Qp);
      %[v,rho,eta]=wtcgls(F,WW,d,np,0,1);
      [v,rho,eta]=wtcgls(F,WW,d,np,stpc,step);
      if (f==25) 
        save temp.mat F WW d v 
	%keyboard
      end
      

   elseif (method==4)
      %Qp=1*ones(np,1);
      %Qp(15)=0.01;
      %Qp(29)=0.01;
      WW=diag(Qp);
      [v,rho,eta,UU,BB,VV] = wlsqr(F,WW,d,np,stpc,reorth,step);
   elseif (method==5)
      [UU,SS,VV]=svd(F);
      SSI=SS.';   %'
      %wmin=1e10*min(diag(SS));
      wmin=1e-1;
      k=0;
      for i=1:min(nh:np)
          if (SS(i,i) > wmin) SSI(i,i)=1/(SS(i,i)+eps);k=k+1;
          else SSI(i,i)=0;
          end,
      end,
      ke(f)=k;
      v=(VV*SSI*UU')*d;  %'
   elseif (method==6)
      [v]=cg0((FH*F+0.01*eye(np)),FH*d);
   elseif (method==7)
      [UU,ss,VV] = csvd(F);
      k_tsvd = gcvnp(UU,ss,d,'tsvd'),%pause;
      if (k_tsvd>=np/3) k_tsvd=np/3;end
      v = tsvd(UU,ss,VV,d,k_tsvd);
   elseif (method==8)
      [UU,ss,VV] = csvd(F);
      %lambda=gcvnp(UU,ss,d,'Tikh'),
      %if(lambda<1e-4) lambda=1e-4;end
      %if(lambda>1e2) lambda=1e2;end
      lambda=1;
      [v,rho,eta] = tikhonov(UU,ss,VV,d,lambda);   
      
   elseif (method==9)
      WW=diag(Qp);      
      [UU,sm,XX] = cgsvd(F,WW);
      k_tgsvd = gcvnp(UU,sm,d,'tsvd'),%pause;
      if (k_tgsvd >=np/4) k_tgsvd=round(np/4);end;
      v = tgsvd(UU,sm,XX,d,k_tgsvd);
   elseif (method==10)
      WW=diag(Qp);      
      [UU,sm,XX] = cgsvd(F,WW);
      %lambda=1;
      lambda=gcvnp(UU,sm,d,'Tikh');
      if(lambda<1e-4) lambda=1e-4;end
      %if (f>2) [v,rho,eta] = tikhonov(UU,sm,XX,d,lambda);
      %else     [v,rho,eta] = tikhonov(UU,sm,XX,d,lambda);
      %end;
      [v,rho,eta] = tikhonov(UU,sm,XX,d,lambda);
    elseif (method==11)
      WW=diag(Qp);
      
      if ((f==2)|(f==10)|(f==30)|(f==(nf-30)))
         [UU,ss,VV] = csvd(F);
         index=index+1;
         subplot(4,1,index);picard(UU,ss,d);
         axis([0 40 1e-10 1e+10]);
      end
      figure(gcf);
      [v,rho,eta]=wtcgls(F,WW,d,np/3,0);  
    elseif (method==12)
      WW=diag(Qp);
      [v,rho,eta]=wtcgls(F,WW,d,np,stpc,step);
      if (f==(80))
         load noise
         noise=E(f,:).'; %' Noise at frequency f
         reguradon(F,d,Vtrue(f,:),WW,noise,10,2,eps1,eps2);
      end
    
   end        
   V=[V; v.'];
   dr=F*v;
   DR=[DR; dr.'];
end
V=seis_shape(V);
size(V)
toc
display('end of loop')
end



if (method==4) figure,plot(w,ke);end,
save condn.mat condn











