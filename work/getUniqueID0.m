function [index,rangey]=getUniqueID0(x,y)
% given two vectors with coordinates find a unique number
% to identify every element and such that the original 
% coordinates can be recovered
% Note: this fails with negative numbers. Use getUniqueId instead.

n=length(x);
if (n ~= length(y)) display('error: vectors x and y need to be the same length');end


miny = min(y);
maxy = max(y);

rangey = maxy - miny + 10;

index =x*rangey + y;

return


