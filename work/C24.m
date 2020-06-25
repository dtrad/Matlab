 t1=-10:.01:0;t1=t1';
g1=-2*exp(2*t1);
t2=0:.01:10;t2=t2';
g2=2*exp(-t2);
t=[t1;t2]; g=[g1;g2];
f=[zeros(size(g1));ones(size(g2))];
c=0.01*conv(f,g);
t=-20:.01:5;t=t';
plot(t,c(1:length(t)))    
