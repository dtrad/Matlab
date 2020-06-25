function [y,ysteps]=newtonparam2(d1fx,d2fx,x0,x,p,tol,itermax);
%input: 
%	
% d1fx : first derivative function  (function for finding roots)
% d2fx : second derivative function (first derivative for finding roots)
% x0 : starting point
% tol: tolerance 
% itermax: maximum iterations
%
% output:
%  y : x cordinate for extreme point of f(x)
%  ysteps :  vector with the history for x as a function of the
%  number of iterations
%
% The same routines can be used to find the minimum or maximum and zeros
% When the zero is required the input arguments are the function
% and  the first derivative.
% When the extremum is required the input arguments are the first
% and 
% second derivative. 
%
% Daniel Trad -- EOSC555 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter=1;
error=0;
y=x0;
ysteps(iter)=y;
x=x(:);
p=p(:);

while( (error>tol | iter==1) & iter < itermax  )
  
    %d1fy=feval(d1fx,y,x,p);
    d1x=feval(d1fx,x+y*p,1);
    d1fy=d1x'*p;
    %d2fy=feval(d2fx,y,x,p);
    d2x=feval(d2fx,x+y*p,2);
    d2fy=p'*d2x*p;
    
%   Loop to move away of x when the second derivative is zero 
%   Hopefully it never goes into this
%   eps is the small number defined in matlab. The factor 100 I
%   used is empprical and it could be defined 
%   as an argument instead to speed up computations. 

    while (abs(d2fy) < eps & abs(d1fy) > tol) 
      if (d1fy > 0 )
	y=y-100*eps
      else 
	y=y+100*eps
      end
      d1x=feval(d1fx,x+y*p,1);
      d1fy=d1x'*p;      
      %d1fy=feval(d1fx,y,x,p);
      %d2fy=feval(d2fx,y,x,p);
      d2x=feval(d2fx,x+y*p,2);
      d2fy=p'*d2x*p;      
    end  
    
    yn=y-d1fy/d2fy;
%   Change to find always a minimum
    %yn=y-d1fy/abs(d2fy);
    %error=abs(feval(d1fx,yn,x,p));
    d1x=feval(d1fx,x+yn*p,1);
    error=abs(d1x'*p);
    iter=iter+1;
    ysteps(iter)=yn;

    y=yn;
      
  

end



