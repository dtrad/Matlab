function [A]=kernel3(nt,dt)
% function [A]=kernel3(nt,dt) 
% it generates a set of cos basis functions with increasing frequency
A=[];
t=0:nt-1;t=t*dt;
for i=1:nt;
  g=cos(2*pi*i*t);
  %plot(g);figure(gcf) % Uncomment to see the basis functions
  A=[A,g(:)];
end

return;


