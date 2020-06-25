
load xxc;h0=(0:49)*dh;
[xxn,dtii]=rnmo3b(xxc(:,1:50),3500,h0,dt,+1); 
xxnn=rnmo3b(xxn(:,1:50),3500,h0,dt,-1,dtii); 

figure,wigb(xxc(1:512,1:50));title('original')
figure,wigb(xxnn(1:512,1:50)); title('nmo and back')