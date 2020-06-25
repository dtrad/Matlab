function [data_out,offset_axis,time_axis,vr,pp_axis,vrorig]=interp_taup(u,h0,h1,dt,pnorm,optionp,vel_min,np,NITER,hyperparameter,vel_max,optionmute,freqint,option_residuals)
% [data_out,offset_axis,time_axis]=interp_taup(u,h0,h1,dt,pnorm,optionp,vel_min,hyperparameter,optionmute,freqint)
% Taup interpolation-  Daniel Trad- 09/08/98
% Interpolation of seismic data from tau-p domain
% Based on Improving Resolution of Radon Operators using 
% a model re-weighted least square procedure 
% Mauricio Sacchi and Tad Ulrych.
% Journal of Seismic Exploration 4, 315-328 (1995)
% Input data
echo on
u=seis_shape(u);
[nt,nh]=size(u);
uorig=u;
if nargin<3|isempty(h1) h1=h0;end
if nargin<4|isempty(dt) dt=0.004;end            	% Sample interval
if nargin<5|isempty(pnorm) pnorm=1;end					% Model Norm to minimize
																	% pnorm=1; -->L1 norm (Huber)
																	% pnorm=2  -->L2 norm (least squares)
																	% pnorm=10; --> Cauchy's norm. 

if nargin<6|isempty(optionp) optionp='linear ';end		% 'linear' slant stack (stack along straight lines), 
																	% 'Yilmaz' stack along hyperboles with time stretching
																	% 'Hampson' stack along paraboles on residual nmo.
                                      
if nargin<7|isempty(vel_min) vel_min=1000; end   	% Minimum velocity                                          
if nargin<8|isempty(np) np=50;end     					% Number of traces in Tau-p domain
if nargin<9|isempty(NITER) NITER=5;end 					% Maximum number of iterations for non-linear problem
if nargin<10|isempty(hyperparameter) hyperparameter=0.1;end  % Hyperparameter for inversion  
if nargin<11|isempty(vel_max) vel_max=5000;end  % Velocity to use in NMO
if nargin<12|isempty(optionmute) optionmute='n';end  % Velocity to use in NMO
if nargin<13|isempty(freqint) freqint=[2 nt/2-10];end  % Velocity to use in NMO

% Settings for inversion
tol=0.01; 			% Tolerance to use in iterations
tolmodel=1e-1;  	% Relatives differences in norm less than this will stop iterations
epsilon=1e-4;   	% Minimum value for v(ii) in L1 (Huber parameter) 

if optionp=='linear ' 
   Power=1;
elseif (optionp=='Yilmaz '|optionp=='Hampson') 
   Power=2;
end;

% Take only part of the data to try interpolation

[HH,WU0,WU1]=parameters_tau_p(h0,h1);

nh0=length(h0);
nh1=length(h1);

h_near=h0(1);
clear data_in;
dts=dt.*2; 		% Only for Yilmaz


if optionp=='linear '
pmax=1/vel_min; 			% Maximum p
p0=0 				 	% Minimum p
elseif optionp=='Hampson'
   max_h=max(h1);
   min_h=min(h1);
   if min_h==0 min_h=h1(2)-h1(1);end
   pmax=1/max_h^2;
   p0=-pmax;
elseif optionp=='Yilmaz '
pmax=1/vel_min;
p0=1/vel_max;
end

% Residual nmo
if optionp=='Hampson'  
   [unmo,dtii]=rnmo3b(u,vel_max,h0,dt,+1);
   u=unmo;clear unmo;
end;

% compute the covariance of the noise (in freq)
%randn('seed',1); % Change seed for noise 
%noise=0.01*randn(nt,nh0);
%u=u+noise;
%clear noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the covariance of the noise (in freq)
randn('seed',2); % Change seed for noise 
noise=0.01*randn(nt,nh0);
NOISE=fft(noise);clear noise;
C=cov(NOISE);clear NOISE;
% In the general case CI=inv(C); but because I assume white noise
CI=diag(1./diag(C));
% Thus CI is a diagonal matrix containing the standar deviation of noise
% W contains the variance of noise
W=diag((diag(CI).^0.5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=normalize(u);

% Offset, velocity axis and Time axis
t=0:dt:((nt-1)*dt);% Time axis.
[alfa,dalfa]=param(optionp,p0,pmax,np);
uini=u; % Save original data
[message]=alias_check(alfa,h1,60)
% Stretching

if optionp=='Yilmaz '  
   [us]=t_2(dt,dts,u,1); %1 means t --> t^2
   dtini=dt;dt=dts;
	u=us(1:length(us(:,1)),1:nh);clear us;
end

% Freq-Time transformation
U=fft(u,nt); % data are zeropad to nt power of 2.
fs=1/dt;  % sample freq.
w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.
UD=duplic(U(1:nt/2,:));
ur=ifft(UD); % Reconstructed data after Fourier Transform 

if optionp=='Yilmaz '
[tmp]=t_2(dtini,dt,ur,-1);dtini,dt,clear ur;
for ii=1:nh;ur(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First reconstruction:(no inversion) 
% a)Forward Transform:t,x--->tau,p 
% b)Backward Transform:tau,p-->t,x
figure,
% I work with half of the FT because the traces are real.
UH=U(1:nt/2,:);   

% h:offset--> row vector
% p:slowness--> column vector

h0=h0(:)';
h1=h1(:)';

% In a general case WU and WV account for irregular geometry in x and p.
% I put the simplest situation.

WV=dalfa*eye(np);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non linear Inversion
% Initial Model
ii=1:np;
V=zeros(nt/2,np);

[V,JPF,J,JD,JP]=forwnl3(UH,V,epsilon,NITER,tolmodel,pnorm,hyperparameter,C,w,h0,nt,Power,WV,WU0,W,np,alfa,freqint);

V=zeropadm(V,nt/2);
VD=duplic(V);
vr=ifft(VD);
   
plotcost(JPF,J,JD,JP)

if optionmute=='y'
   display('Use vr=mute_matrix(vr) (many times) to mute rectangles in vr') 
   display('Then write return and press enter')
   vrorig=vr;
keyboard
V=fft(real(vr));
V=V(1:nt/2,:);
end

[utr]=backward2(V,UH,w,h1,nt,Power,WV,WU0,WU1,W,alfa);

if optionp=='Yilmaz '
    [tmp]=t_2(dtini,dt,utr,-1);dtini,dt,clear utr;
    for ii=1:nh1;utr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
	 [tmp]=t_2(dtini,dt,vr,-1);dtini,dt,clear vr;
    for ii=1:np;vr(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
    if optionmute=='y'
    [tmp]=t_2(dtini,dt,vrorig,-1);dtini,dt,clear vrorig;
	 for ii=1:np;vrorig(:,ii)=paddyad(tmp(:,ii));end;clear tmp;
	 end
end;
ur=uorig(1:nt,:);
utr=utr(1:nt,:);
vr=vr(1:nt,:);

if Power==1 	 				
   		figure;plotwig1(ur,vr,utr,t,h0,h1,alfa,HH)
      elseif Power==2
         if optionp=='Hampson'
            ur=rnmo3b(ur,vel_max,h0,dt,-1,dtii);
            [tempp,dtii]=rnmo3b(utr,vel_max,h1,dt,+1);
            utr=rnmo3b(utr,vel_max,h1,dt,-1,dtii);
         end;
   
         figure;plotwig2(uorig(1:nt,:),vr,utr,t,alfa,h0,h1,HH,option_residuals)
end   

data_out=utr;clear utr;
offset_axis=h1;
time_axis=t;
pp_axis=alfa;

%mytitle2(pnorm);
%=====================================================================    
 


