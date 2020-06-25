function [m]=minphi(A,Wm,Wd,t,beta,m0,tol,maxit,option)
% Minimize the cost function using conjugate gradient
% Solving the system of equations 
% b=Ax
% where 
%		b=G'*Wd'*Wd*d+(Wm'*Wm)*m0
%     A=G'*Wd'*Wd*G + Wm'*Wm
%     x= model
%     m0=ref model
% Daniel Trad - Geop507
[nr nc]=size(A);

if (nargin<6|isempty(m0)) m0=zeros(nc,1);end
if (nargin<9|option==[]), option='in';end
if (nargin<8|maxit==[]), maxit=length(t);end
if (nargin<7|tol==[]), tol=1e-7;end

Ad=Wd*A;
WmtWm=Wm'*Wm;
rhs=Ad'*Wd*t+WmtWm*m0;
lhs=Ad'*Ad+beta*(WmtWm);

if (option=='cg') m=cgs(lhs,rhs,tol,maxit);
else 
   m=inv(lhs)*(rhs);
end;

