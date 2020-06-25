function [t]=dsr(t0,h,v,xm)

t=sqrt(t0.^2/4 + ((xm+h)/v).^2) + sqrt( t0.^2/4 + ((xm-h)/v).^2) ;

