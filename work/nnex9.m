close all
clear all
P=randn(2,16);
plot(P(1,:),P(2,:),'+r');

net=newc([0 1;0 1],8,.1);
w=net.IW{1};
hold on;  h = plot(w(:,1),w(:,2),'ob');
net.trainParam.epochs=40;

net=train(net,P);

w=net.IW{1};
hold on;
delete(h);
h = plot(w(:,1),w(:,2),'ob');
p = [0; 0.2];
a = sim(net,p)