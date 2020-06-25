%Each of the five column vectors  in P defines a 2-element input vectors:

P = [ -0.5 -0.5 +0.3 -0.1;
      -0.5 +0.5 -0.5 +1.0]

%A row vector T defines the vector's target categories.
T = [1 1 0 0];

%Plot these vectors with PLOTPV:
plotpv(P,T);

%NEWP creates a network object and configures it as a perceptron.
net=newp([-1 1; -1 1],1);


plotpv(P,T)

linehandle=plotpc(net.IW{1},net.b{1})

%Train
E=1;
while (sse(E))
    [net,Y,E]=adapt(net,P,T);
    linehandle=plotpc(net.IW{1},net.b{1},linehandle); drawnow;
end;

%New vector
p=[0.7; 1.2];  a = sim(net,p);
plotpv(p,a);
ThePoint=findobj(gca,'type','line');
set(ThePoint,'Color','red');

hold on;
plotpv(P,T)
plotpc(net.IW{1},net.b{1});
