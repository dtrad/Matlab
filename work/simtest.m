
%P3=1:0.05:1.45;
P3=P3(:);
P3=P3(1:10,1);
figure(7),plot(P(:),T(:),'o',P3,sim(net,P3),'+');
figure(gcf)


