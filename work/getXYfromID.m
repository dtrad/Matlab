function [x,y]=getXYfromID(index,rangey,minx,miny)
% undo lexicographic ordering in two dimensions
% given a unique ID number and the range of the last one
% return the two original coordinates  
% Daniel Trad
x = zeros(size(index));
y = zeros(size(index));

x = floor((index - 1)/rangey)+minx;
y = index - 1 - (x-minx)*rangey + miny;

return