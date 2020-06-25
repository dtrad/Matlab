function [d]=readsudata(name,nt,nx)
% read seismic unix data.
% Usage  [d]=readsudata(name,nt,nx)
% 
fid=fopen(name,'r')
d=fread(fid,inf,'float32');
d=reshape(d,nt,nx);
fclose(fid);

return