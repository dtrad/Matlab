function [y,ysteps]=newtonparam3(fx,x0,x,p,b,A,AA,lambda,tol,itermax);
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
%display('inside line search');
while( (error>tol | iter==1) & iter < itermax  )
  
    d1x=feval(fx,x+y*p,b,A,AA,lambda,1);
    d1fy=d1x'*p;
    d2x=feval(fx,x+y*p,b,A,AA,lambda,2);
    d2fy=p'*d2x*p;
    
%   Loop to move away of x when the second derivative is zero 
%   Hopefully it never goes into this
%   eps is the small number defined in matlab. The factor 100 I
%   used is empprical and it could be defined 
%   as an argument instead to speed up computations. 

    while (abs(d2fy) < eps & abs(d1fy) > tol) 
      display('Inside slow loop....');
      if (d1fy > 0 )
	y=y-100*eps
      else 
	y=y+100*eps
      end
      d1x=feval(fx,x+y*p,b,A,AA,lambda,1);
      d1fy=d1x'*p;      
      d2x=feval(fx,x+y*p,b,A,AA,lambda,2);
      d2fy=p'*d2x*p;      
    end  
    
%   Change to find always a minimum
    yn=y-d1fy/abs(d2fy);
    d1x=feval(fx,x+yn*p,b,A,AA,lambda,1);
    error=abs(d1x'*p);
    iter=iter+1;
    ysteps(iter)=yn;
    y=yn;
    
end
iter;
%display('exiting line search');







