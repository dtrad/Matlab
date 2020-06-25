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
global w dte p h t np nh nt Power WV WU W dh h_near ff tor dp dt alfa


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
optionadj='y'  % Reconstruction through the adjunt

% 'linear' slant stack (stack along straight lines), 
% 'Yilmaz' stack along hyperboles with time stretching
% 'Hampson' stack along paraboles on residual nmo.
 

%optionp='linear';
optionp='Yilmaz ';
%optionp='Hampson';
optionmute='n';

if optionp=='linear ' 
   Power=1;
elseif (optionp=='Yilmaz '|optionp=='Hampson') 
   Power=2;
end;

if (optionc==3|optionc==0) pnorm=2;end, 

load c:\daniel\taup521\output.mat
u=xxc(1:512,11:50);
clear xxc;

%Field Parameters 
min_depth=600;
max_depth=3600;
vel_min=2000;
vel_max=4000;

[nto nh]=size(u); % Original data set, nh:Number of traces 

nt=512; % Lenght of time series to work with FFT
dt=0.004; % sample interval
if optionp=='Yilmaz ' dts=dt.*2;end % Only for Yilmaz
dte=1; % Additional factor for Basis functions
t0=0;  % Initial time
dh=50; % Offset
h_near=500; % Closest trace
axis1=[0 (nt*dt) h_near (h_near+dh*nh)];
np=40;     % Number of traces in Tau-p domain
           
         
%pmax=sin(atan((h_near+dh*nh)/min_depth))/vel_1 % Maximum p
if optionp=='linear '
pmax=1/vel_min; % Maximum p
p0=sin(atan((h_near+dh*1)/max_depth))/vel_max; % Minimum p
elseif optionp=='Hampson'
pmax=1/vel_min;
p0=1/vel_max;
elseif optionp=='Yilmaz '
pmax=1/vel_min;
p0=1/vel_max;
end

% Settings for mute 
if optionmute=='y'
tim=0.1;
tfm=0.59;
u=mute(u,tim,tfm);
end

% Residual nmo
if optionp=='Hampson'  
   unmo=rnmo(u,nto);
   u=unmo;
end;

% Settings for inversion
initer=1;  
enditer=1; 
CHIFACT=1;
[MM NN]=size(u);targetphi=MM*NN*CHIFACT;
NITER=5; % Maximum number of iterations for non-linear problem

% Values to use in sigmap
for kk=initer:enditer;sigmap(kk)=(10^(kk-3));end

tol=0.00001; % Tolerance to use in iterations
tolmodel=1e-2;  % Relatives differences in norm less than this will stop iterations
epsilon=1e-6;   % Minimum value for v(ii) in L1 (Huber parameter) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Noise
randn('seed',200); % Change seed for noise 
noise=0.0001*randn(nto,nh); % Noise to add to the shyntetic data
u=u+noise;  % data--->data+noise

% compute the covariance of the noise (in freq)
randn('seed',2); % Change seed for noise 
noise=0.0001*randn(nto,nh);
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
t=t0:dt:((nt-1)*dt);t2=t.^2; % Time axis.
h =h_near:dh:((nh-1)*dh+h_near); % Offset axis

dp = (pmax-p0)/(np-1);
p = p0:dp:pmax; 

if optionp=='Yilmaz '  
   alfa=(p.^2).';
   [us]=t_2(dt,dts,u,1); %1 means t --> t^2
   uini=u;dtini=dt;dt=dts;
   u=us(1:nt,1:nh);clear us;
end

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
if optionp=='Yilmaz '
[tmp]=t_2(dtini,dt,ur,-1);dtini,dt,
for ii=1:nh;urr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
end;

% Plots- 
plotdat0(uini,U,urr,axis1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First reconstruction:(no inversion) 
% a)Forward Transform:t,x--->tau,p 
% b)Backward Transform:tau,p-->t,x
figure,
% I work with half of the FT because the traces are real.
UH=U(1:nt/2,:);   

% h:offset--> row vector
% p:slowness--> column vector
h=h(:)';%[MM NN]=size(h); if MM>NN h=h.';end;
p=p(:) ;%[MM NN]=size(p); if MM<NN p=p.';end;
% In a general case WU and WV account for irregular geometry in x and p.
% I put the simplest situation.
WU=dh*eye(nh);
WV=dp*eye(np);

if optionadj=='y'
% Forward Transform: v=Lu or v=FWU.u

[V]=forward(UH);
V=V.';
VD=duplic(V);
vr=ifft(VD);

% Backward Transform: u~=L*v or u~=F*WV.v
% u~ are the recovered data, which are not equal to u
[utr,JDFB]=backward(V,UH);

if optionp=='Yilmaz '
[tmp]=t_2(dtini,dt,utr,-1);dtini,dt,
for ii=1:nh;utr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
[tmp]=t_2(dtini,dt,vr,-1);dtini,dt,
for ii=1:nh;vr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
end;

temp=energy(utr,urr); % energy of utr equal to energy of data 
plotdata(urr,vr,temp,axis1);
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse problem

% If pnorm=2 --> linear problem--> no iteration.
% pnorm=1 (L1) and pnorm=3 (Cauchy norm) --> nonlinear problem--> iterations

% Linear Inversion
if (optionc==2|optionc==3) % Options with initial model (2) or simplest LS (3)
for kk=initer:enditer; 
		figure,   
      Qpd=1/sigmap(kk).^2;
		Qp=Qpd*eye(np,np);      
      
      [V]=invforw(UH,Qp,CI); % Simplest inversion (Linear)
		
		V=V.';
		VD=duplic(V);
		vr=ifft(VD);
      
      [utr,JD(kk)]=backward(V,UH);
      
      if optionp=='Yilmaz '
		[tmp]=t_2(dtini,dt,utr,-1);dtini,dt,
		for ii=1:nh;utr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
		[tmp]=t_2(dtini,dt,vr,-1);dtini,dt,
		for ii=1:nh;vr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
		end;

      
      JP(kk)=sum(dp*(sum(abs(V(1:60,:).^2))));
      J(kk)=JP(kk)+JD(kk);
      
      if Power==1 
         plotdata(urr,vr,utr,axis1)
      elseif Power==2
         plotdat2(urr,vr,utr,axis1)
      end   
      
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

for kk=initer:enditer;
      kk
      figure,
     % Initial Model
      ii=1:np;
      tol(kk)=epsilon;
      if (optionc==1) V=zeros(size(V));end
      [V,JPF]=ForwNL(UH,V,tol(kk),NITER,tolmodel,pnorm,sigmap(kk),CI);
		V=V.';
		VD=duplic(V);
      vr=ifft(VD);
      [utr,JD(kk)]=backward(V,UH);
      if optionp=='Yilmaz '
			[tmp]=t_2(dtini,dt,utr,-1);dtini,dt,
			for ii=1:nh;utr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
			[tmp]=t_2(dtini,dt,vr,-1);dtini,dt,
			for ii=1:nh;vr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
		end;

      if pnorm==10
         JP(kk)=sigmap(kk).*sum(sum(dp.*log(1+(abs(V(1:60,:))./sigmap(kk)).^2)));
      elseif pnorm==1
         JP(kk)=sigmap(kk).*sum(sum(dp.*(abs(V(1:60,:))./sigmap(kk)).^pnorm));
      else
         JP(kk)=sigmap(kk).^2.*sum(sum(dp.*(abs(V(1:60,:))./sigmap(kk)).^pnorm));
      end   
      J(kk)=JP(kk)+JD(kk);
      if Power==1 
         plotdata(urr,vr,utr,axis1)
      elseif Power==2
         plotdat2(urr,vr,utr,axis1)
      end   
		mytitle(sigmap(kk),pnorm);
end;
     	
plotnorm(sigmap,J,JP,JD);
mytitle2(pnorm);
%=====================================================================    
end;
 
sigmapfin = input('final sigmap?????');

%sigmapfin=findbeta(sigmap,JD,targetphi); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non linear Inversion for the final sigmap
if (optionc~=3)
   initer=1;enditer=1;sigmap(1)=sigmapfin;
   for kk=initer:enditer;
      kk
      figure,

     % Initial Model
      ii=1:np;
      tol(kk)=epsilon;
      if (optionc==1) V=zeros(size(V));end
      [V,JPF]=ForwNL(UH,V,tol(kk),NITER,tolmodel,pnorm,sigmap(kk),CI);
		V=V.';
		VD=duplic(V);
      vr=ifft(VD);
      [utr,JD(kk)]=backward(V,UH);
      if optionp=='Yilmaz '
			[tmp]=t_2(dtini,dt,utr,-1);dtini,dt,
			for ii=1:nh;utr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
			[tmp]=t_2(dtini,dt,vr,-1);dtini,dt,
			for ii=1:nh;vr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
		end;

      if pnorm==10
         JP(kk)=sigmap(kk).*sum(sum(dp.*log(1+(abs(V(1:60,:))./sigmap(kk)).^2)));
      elseif pnorm==1
         JP(kk)=sigmap(kk).*sum(sum(dp.*(abs(V(1:60,:))./sigmap(kk)).^pnorm));
      else
         JP(kk)=sigmap(kk).^2.*sum(sum(dp.*(abs(V(1:60,:))./sigmap(kk)).^pnorm));
      end   
      J(kk)=JP(kk)+JD(kk);
      if Power==1 
         plotdata(urr,vr,utr,axis1)
      elseif Power==2
         plotdat2(urr,vr,utr,axis1)
      end   
		mytitle(sigmap(kk),pnorm);
	end;
 
   figure,
   sumiter=sum(JPF(1:60,:));plot(sumiter,'o');
   title('SUM of JP from 1 to 60 Hz for each iteracion');
   xlabel('iteracion');ylabel('sum JP'); 
end;    

% End of Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%