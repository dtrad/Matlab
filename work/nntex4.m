time = 1:0.01:2.5;
X = sin(sin(time).*time*10);
P = con2seq(X);
T = con2seq(2*[0 X(1:(end-1))] + X);
plot(time,cat(2,P{:}),time,cat(2,T{:}),'--')
title('Input and Target Signals')
xlabel('Time')
ylabel('Input \_\_\_  Target \_ \_')

net = newlin([-3 3],1,[0 1],0.1);

[net,Y,E,Pf]=adapt(net,P,T);
plot(time,cat(2,Y{:}),'b',time,cat(2,T{:}),'r',time,cat(2,E{:}),'g',[1 2.5],[0 0],'k'))

