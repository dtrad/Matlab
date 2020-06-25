P = [1.0   2.0   3.0;
     4.0   5.0   6.0];

T = [0.5 1.0 -1.0];
maxlr=maxlinlr(P,'bias');
net=newlin([0 10;0 10],1,[0],maxlr);
net.trainParam.show = 50; 
net.trainParam.epochs = 500; 
net.trainParam.goal = 0.001; 
[net,tr]=train(net,P,T);
p = [1.0; 4];
a = sim(net,p)