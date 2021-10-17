function [x]=karmarkar_interface(b,A,iter_end,PI_t,DI_t,DG_t,gamma, ...
				 delta,tol,maxit)
% Interface to Log Barrier algorithm karmarkar.m 
%
%   [x]=karmarkar_interface(b,A,iter_end,PI_t,DI_t,DG_t,gamma,delta)
%
% It redefines the non positive model space to a positive space   
% Parameters are the same as karmarkar.m
%
% It finds the solution to the following perturbed LP problem 
% min ||m||_1 + 1/2 || gamma x ||^2 + 1/2 || p ||^2  
% subject to Ax + delta p = b >= 0
% PI_t,DI_t,DG_t are the primary, dual and duality gap tolerance
% gamma and delta are hyperparameters. 
% delta = 1 solves the minimum residual problem
% gamma large aproaches || x ||_2 
% Reference: Chen S. S., Donoho D. L., Saunders M. A. 1998
% Atomic Decomposition by Basis Pursuit, SIAM Journal on Scientific 
% Computing. Volume 20, Number 1, pp. 33-61
%
% Daniel Trad - UBC 2001. 

if (nargin <4) PI_t=1.e-3; end
if (nargin <5) DI_t=1.e-3; end
if (nargin <6) DG_t=1.e-3; end
if (nargin <7) gamma=1.e-4; end
if (nargin <8) delta=1.e-4; end
if (nargin <9) tol=1.e-6; end
if (nargin <10) maxit=100; end

AA=[A -A];
AAS=sparse(AA);
[mm]=karmarkar(b,AAS,iter_end,PI_t,DI_t,DG_t,gamma,delta,tol,maxit);
n=length(mm);
x=mm(1:n/2)-mm(n/2+1:end);





