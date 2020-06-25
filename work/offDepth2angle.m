function [angle]=offDepth2angle(depth,rd,offset)
% give angle for each offset and depth.
% todo: convert to matrix operation
angle = atan((offset/(2*depth-rd)))*180/pi;
return