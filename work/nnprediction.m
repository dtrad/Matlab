clear
clean

ni=20;
no=20;

S2=no;
S1=ni/2;

x=0:ni*no-1;

%P=reshape(x,ni,no);
%P=P'

P=randn(no,ni);

T=(P.^3) - 1* P.^2 + 1*rand(no,1)*ones(1,ni);


 
 net = newff(minmax(P),[S1 S2],{'tansig' 'purelin'});
 
 Y = sim(net,P);
 
 figure(1);
 plot(P',T','+',P',Y','o')

 figure(3);subplot(222);imagesc(Y);
 
 net.trainParam.epochs = 400;
 net = train(net,P,T);

 Y = sim(net,P);


 figure(2);
 plot(P',T','+',P',Y','o')

 figure(3)
 subplot(221);imagesc(T);
 subplot(223);imagesc(Y);
 