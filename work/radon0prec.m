function [v,dr,p]=radon0prec(d,h,np,vmin,vmax,dt,method,vtrue,eps2,step,rmethod)

%  [v,dr,p]=radon0(d,h,np,vmin,vmax,dt,method,vtrue,eps2,step,rmethod)
% Input 
%      d data
%      h offset 
%      np number of p traces
%      vmin min velocity 
%      vmax maximmum velocity
%      dt time interval
%      method  (default Weighted Conjugate Least Squares
%      vtrue   model is available or aproximated model
%      eps2  data variance 
%      step  scale factor for WTCGLS algorithm
%      rtmetdho  Linear Radon Transform
% output
%        v Radon model
%        dr recovered data
%        p  Radon axis
%
% Test different algorithms for computing the Radon transforms,
% in particular using prior information to 
% compute the model covariance and increase resolution in p
% Given the data D computes the Parabolic Radon transform
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
% vtrue is a prior estimation of the model
% eps2 is the variance of the model and acts as an hyperparameter
% step is 1 for perfect data, and between 0.5-0.9 for noise data
% rmethod is one 'LRT' or 'PRT'
%  
% Example
% The following is an example to process the data in data.dat
% Example
% >>load dandata.dat   % dandata.dat contains a matrix x, offset h,time t
% >>[nt,nh]=size(x);
% >>h=0:nh-1;h=h*10;
% >>[v,dr,p]=radon0(x,h,40,1000,3000,0.004,3,[],1e-3,.9,'PRT');
% 
% where 
%      x data
%      h offset 
%      40 number of p traces
%      100 min velocity 
%      3000 maximmum velocity
%      0.004  dt
%      3 method (Weighted Conjugate Least Squares
%      []  vtrue model is available or aproximated model
%      1e-3  data variance 
%      .9    step in WTCGLS algorithm
%      'LRT' Linear Radon Transform
%output
%        v Radon model
%        dr recovered data
%        p  Radon axis
%
% To improve resolution call again the program but this time use
% the previous output as aproximate model
% Example: After first call use
%
% >>vtrue=v;
% >>[v,dr,p]=radon0(x,h,40,1000,3000,0.004,3,vtrue,1e-3,.9,'PRT');
% Daniel Trad-- 6-04-98
global stpc reorth eps2 step;

switch method
 case 1
 display('Cholesky factorization with precomputation');
 case 2
 display('Conjugate gradients Least Squares (niter gives reg)');
 case 3
 display('Weighted Conjugate gradients Least Squares')
 case 4
    display('wlsqr (Weigthed BD) with GCV.')
 case 5
    display('Truncated SVD')
 case 6
    display('Conjugate gradients Least Squares')
 case 7
    display('Truncated SVD and GCV')
 case 8
    display('Tikhonov with SVD with GCV.')
 case 9
    display('Truncated GSVD with GCV.')
 case 12
    display('Tikhonov with GSVD with GCV.')
 case 11
    display('Weighted Conjugate gradients Least Squares')
    display('with eigenvalue plot');
 case 12
    display('Weighted Conjugate gradients Least Squares')
    display ('with regularization demo')
  otherwise, 
    display(' ')
 end
d=normalize(d);
d=seis_shape(d);
[nt nh]=size(d);
% Noise
e=0.01*(randn(size(d)));
E=fft(e);ef=E(nt/2-10,:).';
eps1=real(ef'*ef)/nh
%

d0=d;
save noise.mat d0 e E;
d=d+e;



[nt,nh]=size(d);
stpc=0;   %GCV
reorth=1; %reorthogonalization for lsqr

%nt=256;
if nargin<2|isempty(h) h=0:nh-1;end
if nargin<3|isempty(np) np=nh;end
if nargin<6|isempty(dt) dt=0.004;end 
if nargin<7|isempty(method) method=2;end 
if nargin<8|isempty(vtrue) vtrue=zeros(nt,np);end
if nargin<9|isempty(eps2) eps2=1e-3;end
if nargin<10|isempty(step) step=.9;end
if nargin<11|isempty(method) rmethod='PRT';end 

nargin
vmax
dt
method
max(max(vtrue))
eps2
step
rmethod

if (length(vmin)>1) 
  p=vmin;
else
  p=radonaxis2(np,vmin,vmax,h,rmethod);
end

%load paxis.mat
D=fft(d,nt); % data are zeropad to nt power of 2.
fs=1/dt;  % sample freq.
w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.
DH=D(1:nt/2,:);   

[V,DR]=radonprec(DH,w,h,p,method,vtrue,rmethod);
V=zeropadm(V,nt/2);
DR=zeropadm(DR,nt/2);
VD=duplic(V);
DD=duplic(DR);
v=ifft(VD);
dr=ifft(DD);
t=0:nt-1;t=t*dt;
%figure,wigb(d,1,h,t);
%figure,wigb(dr,1,h,t );   
%figure,wigb(v,1,p,t);

clear global stpc reorth eps2 step;











