function [k]=xkmap(M,N,fs)
% This function computes the |k| for the 2D FT
% Input
%   M Number of rows.(dimension in kx direction)
%   N Number of columns (dimension in ky direction)
%   fs sample frequency.
% Output 
%   Matrix |k|=sqrt(kx^2+Ky^2) (size M x N)
% Taken from R. J. Blakely. pag. 396.

 
for l=1:M;
  for m=1:N;

   dkx=2*pi*(fs/M);
   dky=2*pi*(fs/N);
   nyqx=M/2+1;
   nyqy=N/2+1;

   if (l <= nyqx) 
     kx=(l-1)*dkx;
   else
     kx=(l-M-1)*dkx;
   end
    
   if (m <= nyqy) 
     ky=(m-1)*dky;
   else
     ky=(m-N-1)*dky;
   end
    

   k(l,m)=(kx^2+ky^2)^(1/2);
  end,
end,