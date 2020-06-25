function [kz,kx,w,aa,alfa]=obliquity(dt,dh,NF,NH,vel)
% [kz,kx,w]=obliquity(dt,dh,NF,NH,vel)
% Computes the kz (z component of the wavenumber)
% Daniel Trad- UBC.

w=frequency(dt,NF);
kx=frequency(dh,NH);
ww=w(:)*ones(1,NH);
kxx=ones(NF,1)*kx(NH:-1:1).';
kz=real(sign(ww).*((ww/vel).^2-kxx.^2).^0.5);
alfa=asin(dt/dh*vel);
kzm=ww./vel*cos(alfa);
aa=(abs(kz)>1.01*abs(kzm));
alfa=asin(abs(kxx)./(eps+abs(ww))*vel);
%max(max(abs(imag(alfa(2:256,:)))))