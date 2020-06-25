function [x,y]=getXYfromID0(index,rangey)
% undo lexicographic ordering in two dimensions
% given a unique ID number and the range of the last one
% return the two original coordinates 
% Note: this fails with negative numbers. Use getXYfromID instead.

% Daniel Trad
x = zeros(size(index));
y = zeros(size(index));

x = fix(index /rangey);
y = index - x*rangey;

return