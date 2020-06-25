function [p,a,b]=read_real4(name)
fid=fopen(name,'r');
p=fread(fid,inf,'float32');
fclose(fid);
