function [y]=rotafx(x,ang)
 figure,plot(x);figure(gcf)
 g=imag(hilbert(x));
 figure,plot(g);figure(gcf)
 y=x.*cos(ang./180*pi)-g.*sin(ang./180*pi);
 figure,plot(y);figure(gcf)
 