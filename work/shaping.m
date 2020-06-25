function [y,f]=shaping(epsilon,lf)
%Shaping filter
% epsilon regularization
% lf length of the filter to be used
% y: freq domain filtering
% y2: time domain filtering


lx=128;
ii=1:lx;
x=zeros(1,lx);
x(64)=1;
w=ricker(25,0.004);
w=minphase(w);
w=w(:).';
x=x(:).';
nw=length(w);
w=conv(w,x);w=w(nw/2:lx+nw/2-1);

X=fft(x);
W=fft(w);


%f=real(ifft(conj(fft(w)).*fft(x)./(epsilon+conj(fft(w)).*fft(w))));
%f=real(ifft(fft(x)./(epsilon+abs(fft(w)))));
F=(conj(W).*X./(max(epsilon,(conj(W).*W))));

figure
plot(real(conj(W).*(W)))


Y=W.*F;

Y=duplic(Y(1:lx/2));
y=real(ifft(Y));

F=duplic(F(1:lx/2));
f=real(ifft(F));



y2=conv(w,f(1:lf));
y2=y2(1:lx);

figure
subplot(311);plot(ii,x,ii,w,ii,f)
subplot(312);plot(ii,y)
subplot(313);plot(ii,y2)


