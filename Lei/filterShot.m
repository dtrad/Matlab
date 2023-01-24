function [shotout]=filterShot(shotin,fmax,nt,dt);

shotf=fft(shotin);
%imagesc(abs(shotsf));
df=1/(nt*dt)
shotsw=zeros(size(shotf));
ifmax=floor(fmax/df);
shotsw(1:ifmax,:)=shotf(1:ifmax,:);
shotout=real(ifft(duplic(shotsw(1:nt/2,:))));
figure;
subplot(211);imagesc(shotin);v=caxis();caxis([v(1)/20 v(2)/20]);
subplot(212);imagesc(shotout);v=caxis();caxis([v(1)/20 v(2)/20]);
figure(gcf)
return;