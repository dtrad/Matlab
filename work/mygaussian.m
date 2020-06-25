g,t=mygaussian2(0.50,0.05);
x=gaussian2(t,g,.2);

function [x]=gaussian2(t,g,value);
x=g(find(t==value));
return;
end

function [g,t]=mygaussian2(a,b)
close all
N=1000;
t=0:N;
t=t/N;
g=gaussian(t,a,b,1);
plot(t,g);
return
end




