P = -1:.1:1;
T = [-.9602 -.5770 -.0729  .3771  .6405  .6600  .4609 ...
          .1336 -.2013 -.4344 -.5000 -.3930 -.1647  .0988 ...
          .3072  .3960  .3449  .1816 -.0312 -.2189 -.3201];
plot(P,T,'+');
p = -3:.1:3;
a = radbas(p);
plot(p,a);
plot(p,radbas(p)+radbas(p-1.5)+.05*radbas(p+2),'m-');
eg = 0.02;  % sum-squared error goal
sc = 1;       % spread constant
net=newrb(P,T,eg,sc);
plot(P,T,'+');pause
X=-1:.1:1;
Y=sim(net,X);
hold on;  plot(X,Y);  hold off;