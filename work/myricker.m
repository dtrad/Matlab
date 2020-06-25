      function [w]=ricker(nw,f,dt)

% GENERATES RICKER WAVELET OF LENGTH NW
% WAVELET WILL APPEAR TO BE ABOUT NW/3 POINTS LONG
% function [w]=ricker(nw,f,dt)
% 
      z=1.+(nw-1)/2;
      const=pi*f*dt;
      for ii=1:nw
      a=(ii-z)*const;
      arg=a^2;
      w(ii)=(1-2*arg)*exp(-arg);
      end
      