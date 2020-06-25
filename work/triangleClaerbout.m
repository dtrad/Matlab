function [yy]=triangleClaerbout(nr,nx,xx)
pp=boxcarClaerbout(nr,nx,xx);
qq=boxcarClaerbout(nr,length(pp),pp);
yy=qq(1:nx);
return;