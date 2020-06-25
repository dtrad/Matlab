function [Dout] = gazdag(Din,v,nx,nz,nt,dx,dz,dt,conj) 
%
% Gazdag Migration
% in conj==0 Din(nt,nx) is data and Dout(nz,nx) is image
% in conj==1 Din(nt,nx) is image and Dout(nz,nx) is data
% The dot product is ok.
%

% Common to forw/backward operators
 dkx=2*pi/(nx*dx);
 dw=2*pi/(nt*dt);
 kx=[0:dkx:(nx/2)*dkx -((nx/2)*dkx-dkx):dkx:-dkx];
 w =[0:dw:(nt/2)*dw -((nt/2)*dw-dw):dw:-dw];
 index=find(w == 0);
 w(index)=1E-8;
 [kx w] = meshgrid(kx,w);
 ARG = 1.0-(kx.^2*v^2)./w.^2;
 index=find(ARG <= 0.0);
ARG(index) = 0.0;

% Migration 

if conj==0;

 FFTD=fft2(Din); % Fourier transform data to w-k domain
 IMAGE = zeros(nz,nx); % set up IMAGE matrix
 S  = exp(i*w/v.*sqrt(ARG)*dz);

 for iz=1 : nz
 FFTD = FFTD.*S;
  IMAGE(iz,:) = sum(FFTD,1)/nt;
 end
 Dout=real(ifft(IMAGE,[],2)); 
end;

% Modeling 

if conj==1 ;

  S  = exp(-i*w/v.*sqrt(ARG)*dz);
  IMAGE=fft(Din,[],2);
  FFTD = zeros(nt,nx);


for iz=nz:-1:1
     
 FFTD = FFTD+ones(nt,1)*IMAGE(iz,:);
 FFTD = (FFTD.*S);
end
Dout=real(ifft2(FFTD));

end;

return;