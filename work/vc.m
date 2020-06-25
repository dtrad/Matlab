 function [v]=vc(A)
% This function takes a 2-D array and makes it a vector in loxigonic ordering  
% input:  A - a matrix
% Output: v - a vector 

% E.Haber 1.12.94

 [si,sj]=size(A);
 v=zeros(si*sj,1); 
 for i=1:si,
     v((i-1)*sj+1:i*sj)=A(i,:);
 end;

%=============================================
% Matlab version of
% for j=1:sj,
%   for i=1:si,
%       v((i-1)*sj+j)=A(i,j);
%    end;
%  end;
