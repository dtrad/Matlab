function dout = down_constant_v(din,v,dx,dt);
%
% Migration using the phase shift method
% v is contant. 
%
% din: input data
%   v: velocity (MKS)
%  dx: distance between receivers
%  dt: sampling interval
%
% dout: migrated data  
%
% M.D. Sacchi, Nov 1997, Dept. of Physics, UofA.
%        
% This is part of Dave Michealis Reseach Project
% for GEOPH 326.
%
% sacchi@phys.ualberta.ca
%
%

% 2D fft to data

D = fft2(din) ;

nx = min(size(D));
nt = max(size(D));

% Use the correct fft symmetries and avoid fftshift. This is
% important.

knyq = pi/dx; dkx = 2*pi/(nx*dx); kx = [0:dkx:dkx*(nx/2) -(dkx*(nx/2-1):-dkx:dkx)];
wnyq = pi/dt;  dw = 2*pi/(nt*dt);  w = [0:dw:dw*(nt/2) -(dw*(nt/2-1):-dw:dw)];

i = sqrt(-1);     % The imaginary thing..

for l=1:nt
 for j=1:nx
  kzkz = w(l)*w(l)-v*v*kx(j)*kx(j);    % Predicted squared vertical wavenumber 
   if kzkz > 0
    C(l,j) = exp( i*sign(w(l))*dt*sqrt(kzkz));
     else 
    C(l,j) = 0;
   end
  end  
 end

% Downward Continuation Operator

for l=1:nt
 D = D.*C;                       % Push down the wavefield 
  IMAGE(l,:) =  ones(1,nt)*D;    % imaging condition 
   end

image  = real(ifft(IMAGE,[],2)); 
dout = image; 

