function d = spike(nx,nt,dx,dt,f)
% SPIKE  SPIKE(nx,nt,dx,dt,f) 
%        puts a spike in the center of a grid to
%        the spike is convolved with
%        a Ricker wavelet of central freq. f (in Hz).
%        I use this program to test migration
%        algorithms (to compute the impulse response)
%
%  example:
%
%  d = spike(256,512,25,0.004,40), produces a Ricker wavelet
%  of central freq. f=40Hz in the center of grid of 256 traces
%  by 512 samples. The distance between traces is dx=25m and
%  the sampling interbal dt=0.004sec.
%
%
%  M.D.Sacchi, July 1997, Dept. of Physics, UofA.
%        
%  sacchi@phys.ualberta.ca


wave = ricker(f,dt);  nw = max(size(wave));

d = zeros(nt,nx); 


dx = 20; 
x = 0:dx:dx*(nx-1);


n=nt/2

d(n:n+nw-1,nx/2) = wave(:); 

d = d/max(max(d)); 


