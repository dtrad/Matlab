function [y]=agc3(x,lw)
x=normalize(x);
n=rand(size(x))-0.5;
x=x+n*0.00001;
[nt,nh]=size(x);
y=x;
ii=1:nh;
alfa(ii)=sum(abs(x(1:lw,ii)));
alfam=ones(lw/2,1)*alfa(:).';
y(1:lw/2,ii)=x(1:lw/2,ii).*alfam;
for jj=2:nt-lw
    alfaw(ii)=sum(abs(x(jj:lw+jj-1,ii)));
    y(jj+lw/2-1:jj+lw/2-1,ii)= x(jj+lw/2-1:jj+lw/2-1,ii).*abs(alfa./(alfaw+1e-3));   
end
for ii=1:nh;y(:,ii)=normalize(y(:,ii));end



%figure,
%subplot(211);plot(y);
%subplot(212);plot(alfaw);