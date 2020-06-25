function [test]=checkminimum(fx,x,scale,tol)
% function [test]=checkminimum(fx,x,scale,tol)
% Input: 
%       fx: function
%       x: point to check
%       scale: radius to check around
%       tol: Age of the universe in microseconds (just kidding,
%            tolerance of course!)
% Output:
%       test: if =0 the minimum is all right, if not, you are in
%       trouble my friend. 
%       
% Evaluation of the function fx at random points around the minimum
% Daniel Trad        January 2001
test=0
fz=feval(fx,x)
for i=1:100
  deltax=[scale*randn scale*randn];
  fy=feval(fx,x+deltax);
  if (fy+tol < fz) 
    test=1;
    fy
  end  
    
end

return;