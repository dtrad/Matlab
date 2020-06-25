function [wd,dr]=seewt2(d,L,num1,num2,col,par)
par
qmf = MakeONFilter('Daubechies',par);
wd = FWT2_PO(d,L,qmf);

wd0=wd;
wd=zeros(size(wd0));
num1,num2
wd(1:num1,1:num2)=1;

%wd(1:num,1:num)=wd0(1:num,1:num);;


dr = IWT2_PO(wd,L,qmf);

figure(2),imagesc(wd);colorbar;
figure(3),imagesc(dr);
%figure(4);mesh(dr);
figure(5);plot(dr(1,:)');
%caxis([-col col]);
colorbar;

return;