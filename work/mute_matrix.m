function [xm]=mute(x)
figure,wigb(x)
[ii,jj]=ginput(2);
xm=x;
ii(:)=round(ii(:));
jj(:)=round(jj(:));
xm(jj(1):jj(2),ii(1):ii(2))=0;
figure,wigb(xm);

