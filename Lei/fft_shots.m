function shotsw = fft_shots(nr,ns,nt,dt,fmax)

% nr=384; % number of receivers
% ns=12; % number of shots
% nt=3000; % number of time vector
% dt=0.001; % time interval
df=1/((nt-1)*dt); % frequency interval after DFT 
%fmax=20;  % maximum frequency wished to obtain
shots = readRSF3D('shots.rr',nt-1,nr,ns); % read time domain shots data from rsf file
shotsf=fft(shots); % DFT  data
shotsw=shotsf(1:fmax/df,:,:);  % cut data until maximum frequency


% f=[0:df:fmax]'; % frequency of data from DFT
% fi=[1:fmax]';   % frequency of data after interpolation
% nf=length(fi);
% shotsi=zeros(nf,nr,ns);  % shots data in the frequency domain initialization
% 
% for i=1:nr
%     for j=1:ns
%         
%        shotsi(:,i,j)=interp1(f,shotsw(:,i,j),fi,'pchip');  % interpolation using 'pchip'       
%     end
% end

%imagesc(abs(squeeze(shotsi(10,:,:)))); % check when frequency=10Hz, nr and ns. squeeze is to reduce shotsi(10,:,:) to 2D matrix
%colormap(gray)
end