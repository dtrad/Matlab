function [step]=maxrho2(x,dx)
% Given the equation xn=x+step*dx, this function finds the
% step value (<=1) such that all values in xn are >=0;
% function [step]=maxrho2(x,dx)
% 
m=length(x);
[x1,i]=min(x+dx)
if (x1>=0)  % if all are positive step=1 is OK
  step=1;   
else        % Else find the scale factor to make the minimum = 0 
  step=-x(i)/dx(i);  
end

   


























