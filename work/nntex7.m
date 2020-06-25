close(1)
P = [1 2 3 4 5 6 7 8];
T = [0 1 2 3 2 1 2 1];
plot(P,T,'.','markersize',20)
axis([0 9 -1 4]);figure(gcf)
spread = 0.7; net = newgrnn(P,T,spread);
A = sim(net,P);
plot(P,T,'.','markersize',30);figure(gcf);pause
hold on, plot(P,A,'.','markersize',30,'color',[1 0 0])
p = 3.5;
a = sim(net,p);
plot(P,T,'.','markersize',30);figure(gcf);pause
hold on, plot(p,a,'.','markersize',30,'color',[1 0 0])
P2 = 0:.1:9;
A2 = sim(net,P2);
plot(P2,A2,'linewidth',4,'color',[1 0 0]);figure(gcf);pause
hold on, plot(P,T,'.','markersize',30);figure(gcf);pause
