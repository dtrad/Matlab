function [D] = simple_model(nt,nx,sdev);

% Generation of a very simple model 
% nt: number of samples
% nx: number of traces
% sdev: sdev of the additive Gaussina noise

dt = 4./1000; 
w = ricker(40.,dt); 
nw=max(size(w));
D = zeros(nt,nx);
for i=1:nx
 for j=1:nw
  D(nt/2-floor(nw/2)+j,i)  = w(j);
 end 
end 

% Add noise to the model 
NOISE = sdev * randn(nt,nx);
D = D + NOISE;

return
