function [x2]=slant(a,b)
x=zeros(128,128);
x(a,b)=1;
x2=real(FastSlantStack(x));
subplot(221);imagesc(x);
subplot(222);imagesc(x2);
figure(gcf);
return;
