function[m]=mat(vec,jj)
% [m]=mat(vec,jj)
% This function gets a vector v and returns a matrix with dimensions 
% length(vec)/jj,jj. The numbers are ordered in loxigonic ordering
% This function is the oposite to vc i.e if A is a matrix size ky,kx
% then mat(vc(A),kx)=A

% E.Haber 2.12.94                            

 N=length(vec);
 if size(vec,1)~=1, vec=vec';end;

 ii=N/jj; 
 m=zeros(ii,jj);
 for i=1:ii,
     m(i,:)=vec((i-1)*jj+1:i*jj);
 end;

%=========================================================
% Matlab version of fortran code:
% for i=1:ii,
%  for j=1:jj,
%     m(i,j)=vec((i-1)*jj+j);
%  end;
% end;
