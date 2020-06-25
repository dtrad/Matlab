function zz=kriging(c0,c1,a,x,y,z,xx,yy)
% function zz=krigtest(c0,c1,a,x,y,z,xx,yy)
% Example of ordinary kriging given the variogram
% input
%    c0: nugget effect
%    c0+c1: sill
%    a:  range;
%    x,y,z data points
%    xx,yy unknown points
%    zz    result

x=x(:);y=y(:);z=z(:);
n=length(x);m=length(xx);
C=covariance(x,y,c1,a);
C=[C,ones(n,1)];
C=[C;[ones(1,n) 0]];
CI=inv(C);

D=zeros(n,m);
for i=1:m
    D(:,i)=distance(x,y,xx(i),yy(i),c1,a);
end

D=[D;ones(1,m)];
zz=zeros(m,1);
for i=1:m
    w=CI*D(:,i);
    zz(i)=w(1:n)'*z(:);
end
return

function C=covariance(x,y,c1,a)
% build distance map.
n= length(x);
C=zeros(n,n);
for i=1:n
    for j=1:n
        dist=sqrt( (x(i) - x(j))^2 + (y(i)-y(j))^2 );
        C(i,j)= mycov(dist,c1,a);
    end
end
return

function D=distance(x,y,xx,yy,c1,a)
% build distance map.
n= length(x);
D=zeros(n,1);
for i=1:n
        dist=sqrt( (x(i) - xx)^2 + (y(i)-yy)^2 );
        D(i)= mycov(dist,c1,a);
end
return

function cov = mycov(dist,c1,a)
% covariance function.
cov=c1*exp(-3*abs(dist)/a);
return