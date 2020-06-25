function [index,rangey,minx,miny]=getUniqueID(x,y)
% given two vectors with coordinates find a unique number
% to identify every element and such that the original 
% coordinates can be recovered

n=length(x);
if (n ~= length(y)) display('error: vectors x and y need to be the same length');end

minx = min(x);
miny = min(y);
maxy = max(y);

rangey = maxy - miny + 1;

index =((x-minx)*rangey + (y-miny) + 1);

return


