function [d]=for_taup(m,dt,h,q,N);

%
% FOR_TAUP  [d] = FOR_TAUP(m,dt,h,q,N)
%
%        forward taup liear or parabolic transform
%
% input
%      m(nt,nq) the taup- space 
%      dt       sampling in sec
%      h(nh) offset or position of traces in mts
%      q(nq) ray parameters to retrieve or curvature
%            of the parabola if N=2
%          N =1 Linear tau-p
%            =2 Parabolic tau-p
% output
%      d(nt,nh) the data
%
%      M.D.Sacchi, June 1998, Dept. of Physics, UofA.
%
%      sacchi@phys.ualberta.ca
 
nt= max(size(m));
nh = max(size(h));


M = fft(m,[],1);
D = zeros(nt,nh);
i = sqrt(-1);
for if=2:nt/2;
f = 2.*pi*(if-1)/nt/dt;
L = exp(i*f*(h.^N)'*q);
x = M(if,:)';
y = L * x; 
D(if,:) = y';
D(nt+2-if,:) = conj(y)';
end
D(nt/2+1,:) = zeros(1,nh);
d = real(ifft(D,[],1));






  

