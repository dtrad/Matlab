function plotsvdvectors(F,W,v,d,data,datart,h,p,t)
% Plot the singular vectors for the kernel and for the
% preconditioned kernel
% Input 
%      F kernel
%      W model weights
[n,m]=size(F);ii=1:m;
[U,S,V]=svd(F);


figure(1),
subplot(421)
wigb(data,1,h,t);title('(a)');

subplot(422)
wigb(datart,1,p,t);title('(b)');


subplot(423)
plot(ii,real(v));title('(c)')

subplot(424)
plot(ii,imag(v));title('(d)')

subplot(425)
plot(ii,real(V(:,1:3)));title('(e)')

subplot(426)
plot(ii,imag(V(:,1:3)));title('(f)')


FWWI=F*inv(W);
[U,S,V]=svd(FWWI);

subplot(427)
plot(ii,real(V(:,1:3)));title('(g)')

subplot(428)
plot(ii,imag(V(:,1:3)));title('(h)')




