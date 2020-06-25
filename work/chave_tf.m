function [tf]=chave_tf(ramb1,ramb2,FN)
% Given the output files from Chave ramb1 and ramb2
% tf 
if nargin < 3|isempty(FN) FN=0.25;end
f=ramb1(:,1);
f=f*FN;
zr=ramb1(:,2);
zi=ramb1(:,3);
er=ramb1(:,4);
ei=ramb1(:,5);

[R1,P1]=apres(zr,zi,f);
P1=myunwrap(P1);

f=ramb2(:,1);
f=f*FN;
zr=ramb2(:,2);
zi=ramb2(:,3);
er=ramb2(:,4);
ei=ramb2(:,5);

[R2,P2]=apres(zr,zi,f);
P2=myunwrap(P2);

figure,
subplot(211);loglog(f,R1,'o',f,R2,'+');
subplot(212);semilogx(f,P1,'o',f,P2,'+');VV=axis;axis([VV(1) VV(2) 0 90]);

tf=[f(:) R1(:) R1(:)./10 R2(:) R2(:)./10 P1(:) P1(:)./10 P2(:) P2(:)./10];

