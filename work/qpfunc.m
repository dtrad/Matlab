norm=1;
m=0:100;
eps=1e-7;
eps2=1;
eps1=0.01;
if norm==10
P=-log(1+m.^2./eps2)+12;
Qp=1./(m.^2+eps2)+eps;
subplot(211),plot(m,P)
subplot(212),plot(m,Qp)
elseif norm==1
P=exp(-m.^2./eps2);
Qp=eps2./(2*abs(1+m));
subplot(211),plot(m,P)
subplot(212),plot(m,Qp)
end
