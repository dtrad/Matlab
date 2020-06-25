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

optionc=1;			
pnorm=1;

if (optionc==3|optionc==0) pnorm=2;end, 

load c:\daniel\taup521\output.dat
u=output(26:538,11:50);
clear output1;
%Parameters 
min_depth=600;
max_depth=4000;
vel_1=1000;
vel_m=3300;

[nto nh]=size(u); % Original data set
nt=512; % Lenght of time series
dt=0.004; % sample interval
dte=1; % Additional factor for Basis functions
t0=0;  % Initial time
nh=40; % Number of traces 
dh=50; % Offset
h_near=170; % Closest trace
axis1=[0 2 0 2000];
np=60;     % Number of traces in Tau-p domain
Power=2;   % 1- slant stack (stack along straight lines), 
			  % 2- stack along hyperboles.
           
initern=2;  
enditern=5;           
%pmax=sin(atan((h_near+dh*nh)/min_depth))/vel_1 % Maximum p
if Power==1
pmax=1/vel_1;
p0=sin(atan((h_near+dh*1)/max_depth))/vel_m; % Maximum p
elseif Power==2
%   pmax=0.6/(h_near+dh*nh).^2;
%   p0=0;
pmax=1/vel_1;
p0=sin(atan((h_near+dh*1)/max_depth))/vel_m; % Maximum p
end
%p0=0.0/3400; % Minimum p
tim=0.1;
tfm=0.59;
CHIFACT=1;
[MM NN]=size(u);targetphi=MM*NN*CHIFACT;
NITER=5; % Maximum number of iterations for non-linear problem
sigmap(1)=0.4;
sigmap(2)=0.6;
if pnorm==2 NITER=10;end

tolmodel=1e-2;  % Relatives differences in norm less than this 
					 % will stop iterations
epsilon=1e-6;   % Minimum value for v(ii) in L1 (Huber parameter) 
initerl=2;enditerl=2; % Best sigmap for initial model (2)
if (optionc==3)initerl=1;enditerl=2;end; % Different sigmap to try in (3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%u=mute(u,tim,tfm);

% Noise
randn('seed',200); % Change seed for noise 
noise=0.0*randn(nto,nh); % Noise to add to the shyntetic data

unmo=rnmo(u,nto);u=unmo;
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
t=t0:dt:((nt-1)*dt);t=t.^2; % Time axis.
tol=0.00001; % Tolerance to use in iterarations
h =h_near:dh:((nh-1)*dh+h_near); 

dp = (pmax-p0)/(np-1);
p = p0:dp:pmax;p=p.^2;
      


%p=rnmo2(u);
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
   [utr,vr,V,JD,JP,J,sigmap]=lininv(ur,UH,CI,sigmap,pnorm,C,axis1,initerl,enditerl);
end,

if optionc==2 VLS=V;end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non linear Inversion

if (optionc~=3) 
 	[utr,vr,JD,JPF,J,sigmap]=nlininv(ur,UH,V,CI,epsilon,NITER,tolmodel,pnorm,C,axis1,optionc,initern,enditern);
end;
% With the values on Tikhonov curve I select the sigmap that is closest to the target misfit.
%sigmapfin = input('final sigmap?????');
sigmapfin=1.56e-1;
%sigmapfin=findbeta(sigmap,JD,targetphi); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non linear Inversion for the final sigmap
if (optionc~=3)
   initern=1;enditern=1;sigmap(1)=sigmapfin;
 	[utr,vr,JD,JPF,J,sigmap]=nlininv(ur,UH,V,CI,epsilon,NITER,tolmodel,pnorm,C,axis1,optionc,initern,enditern,sigmap);

   figure,
   sumiter=sum(JPF(1:60,:));plot(sumiter,'o');
   title('SUM of JP from 1 to 60 Hz for each iteracion');
   xlabel('iteracion');ylabel('sum JP'); 
end;    

% End of Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%