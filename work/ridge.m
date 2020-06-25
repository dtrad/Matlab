function [x2]=ridge(a,b)
x=zeros(256,256);
x(a,b)=1;
x2=real(Inv_FastOrthoRidgeletTrans(x));
subplot(211);imagesc(x);
subplot(212);imagesc(x2);
figure(gcf);
return;
