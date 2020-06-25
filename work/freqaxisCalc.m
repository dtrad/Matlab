function [f]=freqaxisCalc(dt,nt)
% function [f]=freqaxisCalc(dt,nt)
% return freq axis with same order as FFT
v1=0:nt/2;
v2=nt/2-1:-1:1;
f=[v1 -v2];
f=f/(dt*nt);
return;
