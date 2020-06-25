%Solving for a Linear Network
%NEWLIND  - Solves for a linear layer.
%SIM   - Simulates a neural network.
%Using the above functions a linear neuron is designed to respond to specific inputs with target outputs.

P = [1.0 -1.2];
T = [0.5 1.0];

w_range = -1:0.1:1; b_range = -1:0.1:1;
ES = errsurf(P,T,w_range,b_range,'purelin');
plotes(w_range,b_range,ES);
figure(gcf);
net = newlind(P,T);
A = sim(net,P);
E = T - A;
SSE = sumsqr(E);

plotes(w_range,b_range,ES);
plotep(net.IW{1,1},net.b{1},SSE);
figure(gcf);
p = -1.2;
a = sim(net,p)
