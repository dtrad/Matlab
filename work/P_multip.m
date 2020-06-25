function [PM]=P_multip(A,B)
% [PM]=P_multip(A,B)
% Multiply 3D matrices
[ns,ng,nt]=size(A);
for ii=1:nt
   PM(:,:,ii)=A(:,:,ii)*B(:,:,ii);
end

