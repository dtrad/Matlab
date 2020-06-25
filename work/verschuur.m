clear;
close all;
dt=0.004;
NTM=32;
load c:\daniel\thesis\shot_gathers.mat; 
xx=xx(1:256,1:1:NTM,1:1:NTM);

[nt,ng,ns]=size(xx);

for ii=1:ns,
   temp(:,:,ii)=zeropadm(xx(:,:,ii),2*nt);
end
xx=temp;
clear temp;
data_sg1=xx(1:2*nt,1:NTM,1);
[nt,ng,ns]=size(xx);

%hannw=hanning(ng);
%hannw=ones(nt,1)*hannw(:).';
%t=0:nt-1;t=t.*dt;t=sqrt(t);tt=t(:)*ones(1,ng);
%for ii=1:ns,
%   xx(:,:,ii)=xx(:,:,ii).*hannw;
%end

%xx(:,17:32,:)=0;

XX=shot_gathers2freq(xx);clear xx;
P=shot_gathers_freq2P(XX);clear XX;

%P=permute(XX(1:256,:,:),[2 3 1]);clear XX;

PM=P_multip(P,P);clear P;
XXR=P2shot_gathers_freq(PM);clear PM;

%XXR=ipermute(PM,[3,1,2]);clear PM;

%w=frequency(dt,nt);w=w(1:nt/2);
%ww=(1-i)*sqrt(w/4/pi);ww=ww(:)*ones(1,ng);
%for ii=1:ns,
%   XXR(:,:,ii)=XXR(:,:,ii).*ww;
%end


xxr=shot_gathers_freq2xt(XXR);clear XXR;

%t=0:nt-1;t=t.*dt;t=t.^(-0.5);t(1:2)=1e-5;tt=t(:)*ones(1,ng);
%for ii=1:ns,
%   xxr(:,:,ii)=xxr(:,:,ii).*tt;
%end

%for ii=1:ns,
%   xxr(:,:,ii)=xxr(:,:,ii)./(hannw+1e-10);
%end

w=ricker(70,0.004);for ii=1:NTM,xxx(:,ii)=deconv_freq_domain(xxr(:,ii,1),w,1e-9);end

for ii=1:4:NTM,figure,wigb(xxr(:,:,ii));end

scale=max(max(data_sg1(150:nt,:)))/max(max(abs(xxx(:,:))));
figure,
wigb((data_sg1+scale*xxx(:,:)));title('data+predicted multiples');
figure,
wigb(data_sg1);title('data')



