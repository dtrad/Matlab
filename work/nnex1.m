% Create a feed-forward backpropagation network.
%	Examples
%
%	  Here is a problem consisting of inputs P and targets T that we would
%	  like to solve with a network.
%
P = [0 1 2 3 4 5 6 7 8 9 10];
T = [0 1 2 3 4 3 2 1 2 3 4];
%
%	  Here a two-layer feed-forward network is created.  The network's
%	  input ranges from [0 to 10].  The first layer has five TANSIG
%	  neurons, the second layer has one PURELIN neuron.  The TRAINLM
%	  network training function is to be used.
%
net = newff([0 10],[5 1],{'tansig' 'purelin'});
%
%	  Here the network is simulated and its output plotted against
%	  the targets.
%
Y = sim(net,P);
plot(P,T,P,Y,'o')
%
%	  Here the network is trained for 50 epochs.  Again the network's
%     output is plotted.
%
net.trainParam.epochs = 50;
net = train(net,P,T);
Y = sim(net,P);
plot(P,T,P,Y,'o')

