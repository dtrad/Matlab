% Final Project- Geop 521B- Daniel Trad- 09/04/98
% Inversion of seismic data from tau-p domain
% Based on Improving Resolution of Radon Operators using 
% a model re-weighted least square procedure 
% Mauricio Sacchi and Tad Ulrych.
% Journal of Seismic Exploration 4, 315-328 (1995)

% Main program: fintau6.m (this code)

% Aditional constructed Functions:

% forward.m: forward p-tau transform
% backward.m is the backward p-tau transform 
% invforw.m is the forward p-tau transform with inversion for least squares
%				(Linear problem)
% forwnl.m is the forward p-tau transform with inversion for 
% 				p and Cauchy norms (Linear problem).
% findtau.m: Given the sigmap-JD values compute the sigmap that yields JD closest to target. 
% duplic.m recovers the full fourier transform from half of it.
% energy.m makes two traces of the same energy. 
% plotrace.m plots a seismic trace.
% plotnorm.m plots the norms as a function of the hyperparameter
% plotdata.m combines different plots with plotrace.m
% plotdat0.m similar for fourier transform
% mytitle.m additional title for plotdata
% mytitle2.m additional title for plotnorm

% Data set: inp.dat x-t domain
% 20 traces sampled at 0.004 seconds 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning environment;
clear;
for jj=1:30,close,end,
clear global
global w dte p h t np nh nt Power WV WU W dh h_near ff tor dp dt


% pnorm=1; -->L1 norm (Huber)
% pnorm=2  -->L2 norm (least squares)
% pnorm=10; --> Cauchy's norm. 
% Optionc:0,1,2,3 
% 0 Least Squares (Linear Inversion), Only for test of the non-linear 
%                                     routine fwdnl.m
% 1 Initial model =0 and then L1 or Cauchy (Non-linear inversion) 
% 2 Initial model =LS and then L1, or Cauchy (Non-linear inversion)
% 3 Simplest LS  -->invforw

optionc=3;			
pnorm=2;

if (optionc==3|optionc==0) pnorm=2;end, 

load c:\daniel\output1.dat
u=output1(:,2:21);

%Parameters 
min_depth=300;
max_depth=800;
vel_1=2000;
vel_m=3200;

[nto nh]=size(u); % Original data set
nt=512; % Lenght of time series
dt=0.004; % sample interval
dte=1; % Additional factor for Basis functions
t0=0;  % Initial time
nh=20; % Number of traces 
dh=50; % Offset
h_near=170; % Closest trace
axis1=[0 1 0 1400];
np=40;     % Number of traces in Tau-p domain
%pmax=sin(atan((h_near+dh*nh)/min_depth))/vel_1 % Maximum p
pmax=1/vel_1;
p0=sin(atan((h_near+dh*1)/max_depth))/vel_m % Maximum p
%p0=0.0/3400; % Minimum p
tim=0.1;
tfm=0.59;
CHIFACT=1;
[MM NN]=size(u);targetphi=MM*NN*CHIFACT;
Power=1;   % 1- slant stack (stack along straight lines), 
			  % 2- stack along hyperboles (not inplemented yet)
NITER=10; % Maximum number of iterations for non-linear problem
sigmap(1)=0.001;
sigmap(2)=0.003;
if pnorm==2 NITER=2;end

tolmodel=1e-2;  % Relatives differences in norm less than this 
					 % will stop iterations
epsilon=1e-6;   % Minimum value for v(ii) in L1 (Huber parameter) 
initialiter=2;enditer=2; % Best sigmap for initial model (2)
if (optionc==3)initialiter=1;enditer=2;end; % Different sigmap to try in (3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=mute(u,tim,tfm);

% Noise
randn('seed',200); % Change seed for noise 
noise=0.0*randn(nto,nh); % Noise to add to the shyntetic data
u=u+noise;  % data--->data+noise
% compute the covariance of the noise (in freq)
randn('seed',2); % Change seed for noise 
noise=0.00001*randn(nto,nh);
NOISE=fft(noise);
C=cov(NOISE);
% In the general case CI=inv(C); but because I assume white noise
CI=diag(1./diag(C));
% Thus CI is a diagonal matrix containing the standar deviation of noise
% W contains the variance of noise
W=diag((diag(CI).^0.5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=u./max(max(abs(u))); % Normalize max amplitude=1

% Offset, velocity axis and Time axis
tor=t0:dt:((nto-1)*dt); % Original time axis
t=t0:dt:((nt-1)*dt); % Time axis.
tol=0.00001; % Tolerance to use in iterarations
h =h_near:dh:((nh-1)*dh+h_near); 
dp = (pmax-p0)/(np-1);
p = p0:dp:pmax;

% Because the number of data is not power of 2, I pad with zeros.
% Freq-Time transformation
U=fft(u,nt); % data are zeropad to nt power of 2.
fs=1/dt;  % sample freq.
w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.

% frequency axis
ff=(0:nt/2-1)*fs/nt; % Positive freq
ff(nt/2+1)=fs/2;% Nyquist
ff((nt/2+2):nt)=ff(nt/2:-1:2); % Negative freq

% Function Duplic.m reconstruc the full Fourier Transform from half of it
UD=duplic(U(1:nt/2,:));
ur=ifft(UD); % Reconstructed data after Fourier Transform 

% Plots- 
plotdat0(u,U,ur,axis1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First reconstruction:(no inversion) 
% a)Forward Transform:t,x--->tau,p 
% b)Backward Transform:tau,p-->t,x
figure,
% I work with half of the FT because the traces are real.
UH=U(1:nt/2,:);   

% h:offset--> row vector
% p:slowness--> column vector
[MM NN]=size(h); if MM>NN h=h.';end;
[MM NN]=size(p); if MM<NN p=p.';end;
% In a general case WU and WV account for irregular geometry in x and p.
% I put the simplest situation.
WU=dh*eye(nh);
WV=dp*eye(np);

% Forward Transform: v=Lu or v=FWU.u
[V]=forward(UH);
V=V.';
VD=duplic(V);
vr=ifft(VD);

% Backward Transform: u~=L*v or u~=F*WV.v
% u~ are the recovered data, which are not equal to u
[utr,JDFB]=backward(V,UH);
temp=energy(utr,ur); % energy of utr equal to energy of data 
plotdata(ur,vr,temp,axis1);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse problem
% compute the covariance of the noise (in freq)
NOISE=fft(noise);
C=cov(NOISE);

% In the general case CI=inv(C); but because I assume white noise
CI=diag(1./diag(C));
% Thus CI is a diagonal matrix containing the standar deviation of noise
% W contains the variance of noise

W=diag((diag(CI).^0.5));

% If pnorm=2 --> linear problem--> no iteration.
% pnorm=1 (L1) and pnorm=3 (Cauchy norm) --> nonlinear problem--> iterations

% Linear Inversion
if (optionc==2|optionc==3) % Options with initial model (2) or simplest LS (3)
   for kk=initialiter:enditer; 
		figure,   
      %sigmap(kk)=(10^kk)*mean(sqrt(diag(C)));
      Qpd=1/sigmap(kk).^2;
		%Qp=Qpd*eye(size(C(1:np,1:np)));
	   Qp=Qpd*eye(np,np);      
      [V]=invforw(UH,Qp,CI); % Simplest inversion (Linear)
		
		V=V.';
		VD=duplic(V);
		vr=ifft(VD);
      
      [utr,JD(kk)]=backward(V,UH);
      
      JP(kk)=sum(dp*(sum(abs(V(1:60,:).^2))));
      J(kk)=JP(kk)+JD(kk);
      
      plotdata(ur,vr,utr,axis1)
      mytitle(sigmap(kk),pnorm);

   end;
   axis
   plotnorm(sigmap,J,JP,JD);
   mytitle2(pnorm);
end,
if optionc==2 VLS=V;end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non linear Inversion
if (optionc~=3) 
   for kk=3:5;
      kk
      figure,
  		sigmap(kk)=(10^(kk-1))*mean(sqrt(diag(C)));
      % Initial Model
      ii=1:np;
      tol(kk)=epsilon;
      if (optionc==1) V=zeros(size(V));end
      [V,JPF]=ForwNL(UH,V,tol(kk),NITER,tolmodel,pnorm,sigmap(kk),CI);
		V=V.';
		VD=duplic(V);
      vr=ifft(VD);
      [utr,JD(kk)]=backward(V,UH);
      if pnorm==10
         JP(kk)=sigmap(kk).*sum(sum(dp.*log(1+(abs(V(1:60,:))./sigmap(kk)).^2)));
      elseif pnorm==1
         JP(kk)=sigmap(kk).*sum(sum(dp.*(abs(V(1:60,:))./sigmap(kk)).^pnorm));
      else
         JP(kk)=sigmap(kk).^2.*sum(sum(dp.*(abs(V(1:60,:))./sigmap(kk)).^pnorm));
      end   
      J(kk)=JP(kk)+JD(kk);
      plotdata(ur,vr,utr,axis1)
		mytitle(sigmap(kk),pnorm);
	end;
     	
   plotnorm(sigmap,J,JP,JD);
   mytitle2(pnorm);
end;
% With the values on Tikhonov curve I select the sigmap that is closest to the target misfit.
sigmapfin=findbeta(sigmap,JD,targetphi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non linear Inversion for the final sigmap
if (optionc~=3) 
      figure,
  		% Initial Model
      ii=1:np;
      tolfin=epsilon;
      if (optionc==1) V=zeros(size(V));end
      if (optionc==2) V=VLS;end;
      [V,JPF]=ForwNL(UH,V,tolfin,NITER,tolmodel,pnorm,sigmapfin,CI);
		V=V.';
		VD=duplic(V);
      vr=ifft(VD);
      [utr,JDfin]=backward(V,UH);
      if pnorm==10
         JPfin=sigmapfin.*sum(sum(dp.*log(1+(abs(V(1:60,:))./sigmapfin).^2)));
      elseif pnorm==1
         JPfin=sigmapfin.*sum(sum(dp.*(abs(V(1:60,:))./sigmapfin).^pnorm));
      else
         JPfin=sigmapfin.^2.*sum(sum(dp.*(abs(V(1:60,:))./sigmapfin).^pnorm));
      end   
      Jfin=JPfin+JDfin;
      %sprintf('The data misfit is %10.1f',JDfin)
      JDfin
      plotdata(ur,vr,utr,axis1)
      mytitle(sigmapfin,pnorm);
      figure,
      sumiter=sum(JPF(1:60,:));plot(sumiter,'o');
      title('SUM of JP from 1 to 60 Hz for each iteracion');
      xlabel('iteracion');ylabel('sum JP'); 
end;    

% End of Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%