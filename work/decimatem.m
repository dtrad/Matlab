function [y]=decimatem(x,R)
% [y]=decimatem(x,R)
% Decimate a matrix column by column
% Daniel Trad- UBC- 22-08-98
x=seis_shape(x);
[NT,NX]=size(x);
y=decimate(x(:,1),R);y=y(:);
for ii=2:NX
   a=decimate(x(:,ii),R);
   y=[y a(:)];
end
   
   