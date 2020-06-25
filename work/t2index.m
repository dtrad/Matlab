function [ii]=t2index(t,dt)
% [ii]=t2index(t,dt)
ii=round(t./dt+1);