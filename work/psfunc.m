function [C]=Condnum(N,W)
%W=1/4;
%lp=0;

%for l=-100:1:100
%   phi(101+l)=sin(2*pi*W*(l-lp))/(pi*(l-lp));
%end
%plot(phi);

alfa=1-cos(2*pi*W);
alfa2=1-cos(2*pi*(1/2-W));
gamma=log(1+2*sqrt(alfa)/(sqrt(2)-sqrt(alfa)));
gamma2=log(1+2*sqrt(alfa2)/(sqrt(2)-sqrt(alfa2)));

Cnum=1-sqrt(pi)*2^(9/4)/sqrt(2-alfa)*sqrt(N)*exp(-gamma*N);
Cden=sqrt(pi)*2^(9/4)/sqrt(2-alfa2)*sqrt(N)*exp(-gamma2*N);
C=Cnum/Cden;