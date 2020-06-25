function [y]=agc(x,lw)
x=normalize(x);
n=rand(size(x))-0.5;
x=x+n*0.00001;
[nt,nh]=size(x);
y=x;
for ii=1:nh
	alfa=sum(x(1:lw,ii).^2);
	y(1:lw/2,ii)=x(1:lw/2,ii)*alfa;
	for jj=2:nt-lw
      alfaw=sum(x(jj:lw+jj-1,ii).^2);
      y(jj+lw/2-1:jj+lw/2-1,ii)= x(jj+lw/2-1:jj+lw/2-1,ii)*abs(alfa/(alfaw+1e-3));   
	end
	y(:,ii)=normalize(y(:,ii));
end


%figure,
%subplot(211);plot(y);
%subplot(212);plot(alfaw);