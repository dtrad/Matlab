% wavelet filter for ground  roll 
% Program for filtering and cleaning MT data contaminated by spikes in time domain.
% and noise due to power line
% calls: grfiltro.m---> subcalls: filtro2gr.m
% 	:fmtfiltr.m---> subcalls: filtro3.m
%
% Daniel Trad - 26-10-96 	
% Copy your data as datos.dat. They must be 5 columns (ex ey hx hy hz).
% Set your parameters in the first part of the code	

%------------------------------------------------
%  Module for setting parameters.
%  Set the parameters here.
%------------------------------------------------  
clear all
close all
        option=1 %Huber (mtfilt)
% 	option=2 %Thomson (mtfilt)
%	filtro='Shrink'
%	filtro='wpdeno'
	filtro='mtfilt'
%	filtro='multi1'
%	Family='Haar',Nfilt=2;
 	Family='Daubechies', Nfilt=16;
%	Family='Symmlet', Nfilt=8;
% 	Family='Coiflet', Nfilt=3;

	FS=1 %sample frequency
	LS=9 %coarser level
	tres=2.5 %treshold
	tipo='MAD';
	Nwindows=1; % Can be greater than size of data.
	ncut=1*1024 % Maximum index to be procesed in data.
	inmin=1;
	inmax=1024*Nwindows; % Max index to be processed.
	FREQ='n'; % Data with noise due to power line.
	GRAPHS='y';  % Graphics of Wavelet coefficient series
	GRAPHW='n';  % Graphics of Wavelet Coefficients
	GRAPHWIND='n';
	GRAPHDATA='n';
	FILT='n';    % Passband in time series.
	ntraces=256;
	ntr0=0;
%-------------------------------------------------
%  Module for loading data
%-------------------------------------------------     
	load /home/dtrad/matlab/grdata2.mat;
	ss=ss(1:1024,ntr0+1:ntr0+ntraces);
	sss=ss(:);	
	ncut=length(sss);inmax=ncut;
	size(sss),pause
	tempdat= sss;%(1:ncut,120); % Cut the index greater than ncut
	
	datos=tempdat;
	clear tempdat;
	MMM=max(size(datos));

	if(inmax>MMM) 
	datos(MMM+1:inmax,:)=datos(MMM-1:-1:2*MMM-inmax,:);
	disp('The data will be wraped around');pause;
	end;
	
	x=detrend(datos(inmin:inmax,1));
	%x=normalize(x);
	clear datos;	
%--------------------------------------------------
%  Module for passband filtering (Fourier)
%--------------------------------------------------

	if(FILT=='y')
	  [bf,af]=butter(2,[0.01 0.9]);
	  x=filtfilt(bf,af,x);
	end;
%--------------------------------------------------
%   Module for setting initial variables
%--------------------------------------------------
	n  = max(size(x));
	D=log(n)/log(2);
	t= (1:n)./FS;	
 	QMF  = MakeONFilter(Family,Nfilt)

%--------------------------------------------------
% Plotting initial data.
%--------------------------------------------------
	if(GRAPHDATA=='y')
	  subplot(111); 
	  plot(t,x)
	  text=sprintf(' Observed ex(t) signal');
	  title(text)
	end;

%--------------------------------------------------
% Options in testing (not working)
%--------------------------------------------------
	%wp=wpanalysis(hx, D, QMF);
	%stree=calcstatree(wp, 'Entropy');
	%[btree,vtree]=bestbasis(stree,D);
	%subplot(111);plotbasistree(btree,D,stree);figure(gcf);pause
	%[wp, btree] = wptour('P',hx,D,QMF,'rambl11');

%	ex=NormNoise(ex,QMF);
%	ey=NormNoise(ey,QMF);
%	hx=NormNoise(hx,QMF);
%	hy=NormNoise(hy,QMF);
%	hz=NormNoise(hz,QMF);
%---------------------------------------------------
%    Module for filtering
%---------------------------------------------------	
if(filtro=='Shrink')
  filtro
  [xh,hxwh]=WaveShr1(x,tipo,LS,QMF);
elseif (filtro=='wpdeno')
  filtro
  [xh,bb,st] = wpdenoi1(x,D,QMF);
elseif (filtro=='mtfilt')
  filtro
  xh = grfiltro(x,LS,QMF,option,tres,GRAPHS,GRAPHW,'x');
elseif (filtro=='multi1')
  filtro
  xh = multimt1(x,LS,QMF);
end
%--------------------------------------------------------

	
subplot(311); 
plot(t,x)
text=sprintf(' x  signal');
title(text);vplot=axis;

subplot(312); 
plot(t,xh)
text=sprintf(' Wavelet Reconstruction from x filtered signal');
title(text);axis(vplot);



%--------------------------------------------------------------------
%    Things not running
%--------------------------------------------------------------------

	%subplot(221),specplot(spectrum(ex,256,[],[],FS),FS)
	%subplot(222),specplot(spectrum(exh,256,[],[],FS),FS)

	%subplot(223),plot(t(3900:4096),exs(3900:4096)),title('EX MT signal');
	%xlabel('time');ylabel('amplitude');
	%axis([3900/FS 4096/FS -60 200])

	%subplot(224),plot(t(3900:4096),exh(3900:4096)),title('EX MT signal');
	%xlabel('time');ylabel('amplitude');
	%axis([3900/FS 4096/FS -600 2000])
%----------------------------------------------------------------------
%   auxiliar shaping
%----------------------------------------------------------------------
x=x(:);
xh=xh(:);
%-----------------------------------------------------------------------
%  substracting the noise
%-----------------------------------------------------------------------
pxh=x-xh;
%------------------------------------------------------------------------
subplot(313); 
plot(t,pxh)
text=sprintf('Removed noise from x(t) signal');
title(text);axis(vplot);
figure(gcf);


xh=reshape(xh,1024,ntraces);
load /home/dtrad/matlab/grdata2.mat;
ss=ss(1:1024,1:256);
figure(3);
subplot(111);wigb(ss(:,1:2:256),10)
sss=xh;
figure(4);	
subplot(111);wigb(sss(:,1:2:256),10);
figure(5);
ssss=ss-sss;
subplot(111);wigb(ssss(:,1:2:256),10);
Ndatos=min(MMM,inmax);
if(GRAPHWIND=='y')		
  for i=1:Nwindows,
    a=512*(i-1)+1;
    b=a+512-1;
    ttex=sprintf('Noise removed:ex(t). Window %d',i);
    MM=max(max(sfd(a:b,1:4)));NN=min(min(sfd(a:b,1:4)));
    subplot(111),plot(t(a:b),sfd(a:b,1)),title(ttex);xlabel('time(sec)')
    figure(gcf),pause
  end,
end,
%       clean signal : sf
%-----------------------------------------------------------------------
	
% Built October 04, 1996 
% Daniel Trad, CONICET.CRICYT. Mendoza. Argentina.

   
    
   

















