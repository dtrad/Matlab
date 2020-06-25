MM=128;NN=64;
I=zeros(MM,NN);
I(50,:)=100;
for ii=1:NN;
   I(60,ii)=100*ii;
   I(70,ii)=100*ii^2;
end

[U,S,V]=svd(I);
II=U*S*V';size(S)
pp=diag(S);
subplot(221);image(I)
subplot(222);image(II)
subplot(223);plot(pp);
pp(1:10)

w=1;
h=0:25:250;nh=length(h);
alfa=0:0.001:0.01;np=length(alfa);
alfa=alfa(:);
v(1:np)=0;v(round(np/2))=1;

WU=eye(nh);
WV=eye(np);

sigmap=1;
% Operators
F=exp(i*w*(alfa*(h.^2)));
FH=F';
L=F*WU;
LH=FH*WV;
for ii=1:np   
Qpd(ii)=(1/2/sigmap^2).*(abs(v(ii)+eps)./sigmap).^(-1);
Qp=diag(Qpd);
end
LL=inv(Qp+L*LH)*L;
