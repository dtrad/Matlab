function plotcostf(lambda)
d=10;
m=-2:0.1:10;
m0=-1;
A=2;
r=d-A*m;
Jd=r.^2;
Jm=lambda*abs(m-m0);
J=Jd+Jm;
plot(m,J,m,Jd,'g.',m,Jm,'r+');figure(gcf);