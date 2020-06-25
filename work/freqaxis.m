function [f]=freqaxis(dt,nt)
% function [f]=freqaxis(dt,nt)
f=(-nt+1)/2:(nt-1)/2;
f=f/(dt*nt);
