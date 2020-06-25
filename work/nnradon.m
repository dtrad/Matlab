fid=fopen('../work/hyperbola','r')
P=fread(fid,inf,'float');
fclose(fid);
P=reshape(P,512,65);
P=P(201:2:300,1:65);
[n1,n2]=size(P);
figure,wigb(P);
figure,simage(P);


 T=zeros(size(P));T(25,35)=1;;
 

 net = newff(ones(n2,1)*[-2 2],n1,{'tansig'});
 Y0 = sim(net,P');
 
 net.trainParam.mu=1e-3;
 net.trainParam.epochs = 10;
 return
 net = train(net,P',T);

 Y = sim(net,P');