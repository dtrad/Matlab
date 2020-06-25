function [dout] = pocs_2d(d,T,dt,f_low,f_high, option, perc_i, perc_f, N, a, tol);
%POCS: 2D (x,t) seismic data regularization using POCS
%
% [dout] = pocs_2d(d,T,dt,f_low,f_high, option, perc_i, perc_f, N, a, tol);
%
%  IN   d:      data (x,t)
%       T:      sampling operator
%       f_low:  min freq in Hz in the data
%       f_high: max freq in Hz inthe data
%       option:    = 1 linear threshold schedulde
%                    2 exponential threshold schedulde
%                    3 data adaptive
%       perc_i: percentage of max amplitude for initial threshold   
%       perc_f: percentage of max amplitude for final threshold (perc_f<<perc_i)   
%       N:      maximum number of iterations 
%       a:      a = 1 for data with high SNR or when one wants to fit the noise 
%               a <1  for noise attenuation (try a=0.2 to 0.4)
%       tol:    tolerance (try 0.01). This is to truncate the number of iterations 
%
%  OUT  dout:  reconstructed data
%
%
%  Copyright (C) 2012, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi and A.Stanton
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%

% Size of cube and fft numbers

 [nt,nx1] = size(d);

  nf = 2^nextpow2(nt);
 nk1 = 2^nextpow2(nx1);
 

% T is sample, S is the reinsertion operator 
 
 S = 1.-T;

% Min and max freq indeces for the reconstruction

 k_low = max(floor(f_low*dt*nf) + 1,2);
 k_high = min(floor(f_high*dt*nf) + 1,nf/2);
 
 Dout = zeros(nf,nx1);
 D = fft(d,nf,1);

 for k = k_low:k_high;
  freq(k) = (k-1)/(nf*dt);
  x = squeeze(D(k,:));
  th = th_schedule(option,x,perc_i,perc_f,N);
  y = x;
   iter = 1;
   E1 = 100;
   
   
   
  while (iter <=N & E1>tol);

   Y = fft(y,nk1);
   A = abs(Y); Angle = angle(Y);
   A_start1(iter,:) = A;
   
   I = find(A<th(iter)); A(I)=0; Y = A.*exp(i*Angle);
   y = ifft(Y); y = y(1:nx1);
   yold = y;
   y = a*x + (1-a)*T.*y + S.*y;
   dif1 = y-yold;
   c1 = sum( (abs(dif1(:))).^2);
   c  = sum( (abs(y(:))).^2);
   E1 = c1/c;
   iter = iter + 1;
   % plot(real(y));figure(gcf);[E1,iter,k]
  end

  iter-1;
   
  Dout(k,:) = y;
  Dout(nf-k+2,:)  = conj(y);
  Amp_s(k,:) = A_start1(1,:);
  Amp_f(k,:) = A;
 end;

 dout = ifft(Dout,[],1);
 dout = dout(1:nt,:);
 
 return

function th = th_schedule(option,x,perc_i,perc_f,N);
% Function used by pocs to define the threshold schedule based 
% on parameter option
%
%  option == 1 --> linear
%  option == 2 --> exponential
%  option == 3 --> use real amplitude to define curve
%

  [nx1] = length(x);
  nk1 = 2^nextpow2(nx1);
  
  X = fft(x,nk1);
  A = abs(X); 
  Amax = max(A(:));
  th_i = perc_i*Amax/100;
  th_f = perc_f*Amax/100;
  k = [1:1:N];

 if option==1;
  th = th_i + (th_f-th_i)*(k-1)/(N-1);
 end;
  
 if option==2;
  b = -log(th_f/th_i);
  th = th_i*exp(-b*(k-1)/(N-1));
 end;

 if option==3 ; 
  Amax = max(A(:));
  avec=sort(reshape(A,nk1,1),'descend');
  I = find(avec>perc_i*Amax/100);
  avec(I)=[];
  I = find(avec<perc_f*Amax/100);
  avec(I)=[];
  th = zeros(N,1);
  th(1) = avec(1);
  th(N) = avec(end);
  for j=2:N-1
  th(j)=avec(ceil((j-1)*length(avec)/(N-1)));
  end
 end;

 return;
