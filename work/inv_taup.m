function [m] = inv_taup(d,dt,h,q,type,N,ilow,ihigh,mu);
% 
% INV_TAUP  [m] = INV_TAUP(d,dt,h,q,type,N,ilow,ihigh)
% 
%        Inverse taup parabolic tau-p transform
% 
% input
%      d(nt,nh) seismic traces
%      dt       sampling in sec
%      h(nh) offset or position of traces in mts
%      q(nq) ray parameters to retrieve or curvature
%            of the parabola if N=2
%      type  =1 Least squares  
%            =2 Sparse inversion via Cauchy criterion
%          N =1 Linear tau-p  
%            =2 Parabolic tau-p
%       ilow =  freq. sample where the inversion starts (if you don't know use 2)
%       ihigh=  freq. sample where the inversion ends   (if you don't know nt/2+1)
%               a typical Ricker wavelet of 30Hz has only 15 non-zero freq. samples 
%               to invert ihigh=15, here is where the comp. saving is done...


% output
%      m(nt,nq) the linear or parabolic tau-p domain 
%
%      M.D.Sacchi, June 1998, Dept. of Physics, UofA.
%
%      sacchi@phys.ualberta.ca
 

nt= max(size(d));
nq = max(size(q));
nh = max(size(h));

D = fft(d,[],1);
M = zeros(nt,nq);
i = sqrt(-1);
alpha=1;
for if=ilow:ihigh
f = 2.*pi*(if-1)/nt/dt;
L = exp(i*f*(h.^N)'*q); tr= trace(L'*L); 
y = D(if,:)';
x = L'*y;

if type == 1
   Q =mu;
   x = inv(L'*L+Q) *L'* y; 
else 
   for iter=1:5
   Q = 0.01*tr*diag(1./(0.0001+x.^2));  
   x = inv(L'*L+Q) *L'* y; 
   end
end


M(if,:) = x';
M(nt+2-if,:) = conj(x)';
end
M(nt/2+1,:) = zeros(1,nq);
m = real(ifft(M,[],1));







  

