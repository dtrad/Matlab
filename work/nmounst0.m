% stretching free nmo project.
% Daniel Trad - Veritas June 2003

close all 
clear all
[d]=readsudata('data1.bin',512,63); 
[d2]=readsudata('data1nmo.bin',512,63);
%d=[d; zeros(size(d))];
%d2=[d2; zeros(size(d2))];
[nt nh]=size(d)
D=fft(d); 
D2=fft(d2); 
D3=zeros(size(D2));
method=4;

switch method
 case 1 % interchange amp and phase
  for i=1:nh
    D3(:,i)=abs(D(:,i)).*exp(sqrt(-1)*angle(D2(:,i))); 
  end
  d3=real(ifft(D3));  
 case 2  % shaping filter
  filter=(conj(D2).*D)./(conj(D2).*D2);
  D3=D2.*filter;
  d3=real(ifft(D3));
  %D3=D3.*conj(filter);
 case 3  % high pass filter
  filter=2*[1;-1];
  
  dtemp=zeros(nt+1,nh);
  for i=1:nh
    dtemp(:,i) = conv(d2(:,i),filter);
    dtemp(:,i) = conv(dtemp(1:nt,i),filter(end:-1:1));
    d3(:,i) = dtemp(1:nt,i); 
  end
  D3=fft(d3);
  
 case 4  % high pass filter
 
  for i=1:nh
    D2H = D2(:,i);
    DH  = D(:,i);
    filter=(conj(D2H).*DH)./(conj(D2H).*D2H);
    D3(:,i) = filter(:).*D2H(:);
    
  end
  d3=real(ifft(D3));  
  
 otherwise
  display 'no method';
end

%for i=1:nh
%  D3(:,i)=sqrt(abs((conj(D(:,i)).*D(:,i)))).*exp(sqrt(-1)*angle(D2(:,i)));  
%end


%D3(:,1:1)=0;
%D3(:,nh-1:nh)=0;
figure(1);
subplot(221);
imagesc(abs(ifft(D)))
colorbar

subplot(222);
imagesc(abs(ifft(D2)));colorbar;

subplot(223);
imagesc(d3)
colorbar

subplot(224);
imagesc(imag(ifft(D3)))
colorbar

figure(2);
subplot(221);
imagesc(abs(D))
colorbar

subplot(222);
imagesc(abs(D2));colorbar;

subplot(223);
imagesc(abs(D3))
colorbar

subplot(224);
imagesc(imag(D3))
colorbar

figure(3),wigb(d2)
figure(4),wigb(d3)
figure(5)
ii=1:nt;
subplot(221);
plot(ii,d(:,nh-15),ii,d2(:,nh-15));
%plot(ii,d(:,nh-10),ii,d2(:,nh-10),ii,d3(:,nh-10));

subplot(222);
plot(ii,abs(D(:,nh-15)),ii,abs(D2(:,nh-15)));


subplot(223);
plot(ii,d(:,nh-15),ii,d3(:,nh-15));
%plot(ii,d(:,nh-10),ii,d2(:,nh-10),ii,d3(:,nh-10);

subplot(224);
plot(ii,abs(D(:,nh-15)),ii,abs(D3(:,nh-15)));
