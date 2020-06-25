function [f,o] = ls_inv_filter(w,NF,Lag,mu);  
%
%LS_INV_FILTER Spiking deconv. This porgram computes an inverse
%        filter, which later is used to deconvolve the
%        seismic trace.
%
%  [f,o] = LS_INV_FILTER(w,NF,Lag,mu)
%
%  IN.
%   Compute the inverse filter of a wavelet
%   NF: the wavelet
%   Lag : the position of the spike in the desired output
%         Lag = 1 (min. phase)
%   mu: Prewhitening in %    
%
%  OUT.
%   f: the filter
%   o: the ouput or convolution of the filer with 
%      the wavelet 
%
%   Example:
%
%   w = [2,1];     % the wavelet
%   [f,o] =  ls_inv_filter(w,20,1,2);  
%   figure(1); plot(f);title('Filter')      % Plot filter and output
%   figure(2); plot(o);title('Wavelet \otimes Filter')
%
%  M.D.Sacchi, July 1998, Dept. of Physics, UofA.
%
%  sacchi@phys.ualberta.ca


NW = max(size(w));    % lenght of the wavelet

NO = NW+NF-1;         % Leght of the output      

[mc,mr]=size(w);
if mc <= mr; w = w'; end;

b = [zeros(1,NO)]';  % Desire output 
b(Lag,1) = 1.;       % Position of the spike

C = convmtx(w,NF);   % Convolution matrix 
 
R = C'*C+mu*eye(NF)/100.;    %  Toeplitz Matrix
rhs = C'*b;                  %  Right hand side vector 
f = inv(R)*rhs;              %  Filter 

if nargout == 2
o = conv(f,w);        %  Actual output
end

return

