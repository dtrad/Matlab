clear
clean

ni=10;
no=10;

S2=no;
S1=5;

x=0:ni*no-1;x=x./20;

reg=1

if (reg)
  
  P=reshape(x,ni,no);
else
  P=randn(no,ni);
end


T=0.2*(P.^3) - 1* P.^2 + 0*rand(no,1)*ones(1,ni);


 
 net = newff(minmax(P),[S1 S2],{'tansig' 'purelin'});
 
 Y0 = sim(net,P);
 
 figure(1);
 plot(P',T','+',P',Y0','o')
 net.trainParam.mu=1e-3;
 net.trainParam.epochs = 1000;
 net = train(net,P,T);

 Y = sim(net,P);


 figure(2);
 subplot(211),plot(P',T','+',P',Y','o')

 figure(3)
 subplot(221);imagesc(P);title('input');
 subplot(222);imagesc(T);title('target');
 subplot(223);imagesc(Y0);title('prediction before training');
 subplot(224);imagesc(Y);title('prediction after training');
 
 
 if (reg)
   x=0:ni*no-1;x=x./20+0.5;
   P2=reshape(x,ni,no);
 else
   P2=randn(no,ni);   
 end
   
 T2=0.2*(P2.^3) - 1* P2.^2 + 0*rand(no,1)*ones(1,ni);
 Y2 = sim(net,P2);
 
 figure(2);
 subplot(212);plot(P2',T2','+',P2',Y2','o')

 figure(5)
 subplot(221);imagesc(P2);title('input');
 subplot(222);imagesc(T2);title('target');
 subplot(224);imagesc(Y2);title('prediction without retraining');

 figure(2);