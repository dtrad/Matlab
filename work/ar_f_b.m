function [D_pred] = ar_f_b(D_in,lf,mu,type);
%
%AR_F_B  : Forward and Backward prediction of a complex
%          TS via autoregressive modelling
%
%
% [D_pred] = ar_f_b(D_in,lf,mu,type)
% 
% IN.
%  D_in: input data in a colunm vector
%  lf  : lenght of the ar process (lenght of the prediction
%        filter)
%  mu  : pre-whitening in %
%  type: 1 Forward, -1 Backward
%       
%  OUT. 
%   D_pred: the predicted TS 
%
%  NOTE: This is also the basis for AR model fitting and
%        and paramtric spectral analysis
%
%  M.D.Sacchi, July 1998, Dept. of Physics, UofA.
%  sacchi@phys.ualberta.ca

        
n = max(size(D_in));
D_pred = zeros(n,1);

if type == -1; D_in=flipud(D_in); end;

% Set convolution matrix and rhs term

C = convmtx(D_in,lf);
d = zeros(n+lf-1,1); d(1:n-1) = D_in(2:n);  

% Compute the filter 

f = inv(C'*C+mu*eye(lf))* C'*d; 
aux = conv(f,D_in);

% Do the  prediction 

D_pred(2:n,1) = aux(1:n-1);

if type == -1; D_pred=flipud(D_pred); end;

return
