function plotdata(ishot)
if (nargin < 1) ishot=1;end
load data.mat;
DD=permute(DATAT,[3 1 2]);
clear DATAT;
nfreq=501;
[nf,nr,ns]=size(DD);
D1=zeros(nfreq,nr);
D1(1:nf,:)=DD(:,:,ishot);
DD1=duplic(D1);
d=(abs(ifft(DD1)));
[n m]=size(d);
figure;imagesc(d);v=axis();caxis([v(1)/10 v(2)/10])
%figure,wigb(d,40,1:m,1:n)
figure(gcf);
end