close(1)
clear all
P = [1 2; 2 2; 1 1]';
Tc = [1 2 3];
plot(P(1,:),P(2,:),'.','markersize',30);
for i=1:3, text(P(1,i)+0.1,P(2,i),sprintf('class %g',Tc(i))), end
axis([0 3 0 3]);figure(gcf);pause
T = ind2vec(Tc)
spread = 1;
net = newpnn(P,T,spread);
A = sim(net,P);
Ac = vec2ind(A);
plot(P(1,:),P(2,:),'.','markersize',30);figure(gcf);pause
for i=1:3, text(P(1,i)+0.1,P(2,i),sprintf('class %g',Ac(i))), end
p = [2; 1.5];
a = sim(net,p);
ac = vec2ind(a);
plot(p(1),p(2),'.','markersize',30)
 text(p(1)+0.1,p(2),sprintf('class %g',ac));figure(gcf);pause
