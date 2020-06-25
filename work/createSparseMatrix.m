function [A]=createSparseMatrix(x,y,z)
% given three vectors of the same length, 
% it creates a sparse matrix A, where row and columms are
% given by x,y and values of A are given by z.

% round to integer all vectors 

x=fix(x);
y=fix(y);

x = x - min(x) + 1;
y = y - min(y) + 1;

A = sparse(x,y,z);

return;



