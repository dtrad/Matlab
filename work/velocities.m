a=ones(100,100);
l1=10;l2=40;l3=60;

 a(1:l1,:)=a(1:l1,:).*1000;
 a(l1+1:l2,:)=a(l1+1:l2,:).*1500;
 a(l2+1:l3,:)=a(l2+1:l3,:).*3500;
 [i,j,v]=find(a);
 pp=[v i j];
 save c:\sunt\vel.dat pp /ascii