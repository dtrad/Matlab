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
method=1;

switch method
 case 1  % high pass filter

    t0  = d(:,floor(nh/2));
    t1  = d(:,nh);
    
    
    
    T0 = fft(t0(125:270));
    T1 = fft(t1(287:337));
    
    nt0 = length(T0);
    nt1 = length(T1);
    
    T0B = [T0(1:floor(nt0/2+1));zeros(nt - nt0,1);T0(floor(nt0/2+2):nt0)];
    T1B = [T1(1:floor(nt1/2+1));zeros(nt - nt1,1);T1(floor(nt1/2+2):nt1)];
    
    t0b = ifft(T0B);
    t1b = ifft(T1B);
    %D2H = D2(:,i);
    %DH  = D(:,i);
    %filter=(conj(D2H).*DH)./(conj(D2H).*D2H);
    %D3(:,i) = filter(:).*D2H(:);

  %d3=real(ifft(D3));  
  
 otherwise
  display 'no method';
end

t0b = t0b(1:nt);
t1b = t1b(1:nt);

figure(1)
ii=1:nt;
subplot(211);
plot(ii,t0,ii,t1,ii,t0b,ii,t1b)
