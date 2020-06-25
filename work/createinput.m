
n=50,m=50;
A=rand(n,m);
x=1:m;
x=x(:);
b=A*x;


save input4.mat A b x

save /home/dtrad/seismic/input4a.txt A /ascii
save /home/dtrad/seismic/input4b.txt b /ascii

!cat /home/dtrad/seismic/input4a.txt  /home/dtrad/seismic/input4b.txt >  /home/dtrad/seismic/input4.txt

