function [output]=energy(input,refer)
% This function makes input of the same energy than refer
% Energy:
ENR=sum(sum(abs(refer)));
ENI=sum(sum(abs(input)));
output=input.*(ENR/ENI);
