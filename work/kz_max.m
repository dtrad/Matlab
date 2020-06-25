function [kzm]=kz_max(dt,dh,NF,NH,vel)
% Given the field parameters returns kz_max to avoid spatial aliasing.
% kz < w/c * cos (inc_angle). This condition is satisfied automatically
% for w_Nyquist, but for lower frequencies it is necesary to check it.
% Daniel Trad- UBC

w=frequency(dt,NF);

kx=abs(frequency(dh,NH));
alfa=asin(dt/dh*vel);


ww=w(:)*ones(1,NH);
kxx=ones(NF,1)*kx(NH:-1:1).';
kz=real(sign(ww).*((ww/vel).^2-kxx.^2).^0.5);
kzm=ww./vel*cos(alfa);
alfa=asin(kxx./ww*vel);

