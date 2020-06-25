function [out]=zeros_2d(in)
% [out]=zeros_2d(in)
% Given the 2D spectrum it makes zeros the 1st and 4th to recover a spectrum 
% of  a real data set dipping towards the right side.
% It maps 3rd from 1st quadrant and 2nd from 4th
% Be careful with DC and Nyquist for 4th quadrant. 
% out(1,NC/2+2:NC) and out(NR/2+1,NC/2+2:NC) are not mapped 
% 
% Daniel Trad, UBC- 30/03/98
[NR NC]=size(in);
out=in;
% 3rd quadrant=1st quadrant
out(NR/2:-1:1,NC/2:-1:2)=0;
out(NR/2+2:NR,NC/2+2:NC)=0;


% 2nd quadrant=4st quadrant
out(2:NR/2,NC/2+2:NC)=conj(in(NR:-1:NR/2+2,NC/2:-1:2));
out(:,NC/2+1)=0;