function [zz]=boxcarClaerbout(nb, nx, xx)
bb=zeros(nx+nb);
ny=nx+nb-1;
bb(1)=xx(1);
for i=2:nx
    bb(i)=bb(i-1)+xx(i); % B(Z)=X(Z)/(1-Z)
end
%vv=nx+1:ny
%bb(vv) = bb(vv-1);
for i=nx+1:ny
    bb(i)=bb(i-1);
end
yy(1:nb)=bb(1:nb);
vv = nb+1:ny;
yy(vv) = bb(vv) - bb(vv-nb); % Y(Z)=B(Z)*(1-Z^nb)
yy=yy/nb;
nh=(nb+1)/2;
zz=yy(nh:nh+nx-1);
return
