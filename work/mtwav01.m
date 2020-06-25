% MT wavelet filter 
% Program for filtering and cleaning MT data contaminated by spikes in time domain.
% and noise due to power line
% calls: mtfiltro.m---> subcalls: filtro2.m
% 	:fmtfiltr.m---> subcalls: filtro3.m
%
% Daniel Trad - 26-10-96 	
% Copy your data as datos.dat. They must be 5 columns (ex ey hx hy hz).
% Set your parameters in the first part of the code	

%------------------------------------------------
%  Module for setting parameters.
%  Set the parameters here.
%------------------------------------------------  
%	option=1 %Huber (mtfilt)
 	option=2 %Thomson (mtfilt)
%	filtro='Shrink'
%	filtro='wpdeno'
	filtro='mtfilt'
%	filtro='multi1'
	Family='Haar',Nfilt=2;
 	Family='Daubechies', Nfilt=4;
%	Family='Symmlet', Nfilt=8;
% 	Family='Coiflet', Nfilt=3;

	FS=0.500 %sample frequency
	LS=4 %coarser level
	tres=3.0 %treshold
	tipo='MAD';
	Nwindows=8; % Can be greater than size of data.
	ncut=8*512 % Maximum index to be procesed in data.

	inmin=1;
	inmax=512*Nwindows; % Max index to be processed.
	FREQ='n'; % Data with noise due to power line.
	GRAPHS='n';  % Graphics of Wavelet coefficient series
	GRAPHW='n';  % Graphics of Wavelet Coefficients
	GRAPHWIND='n';
	GRAPHDATA='n';
	FILT='y';    % Passband in time series.
	GAN='n';

	gex=70;cgex=(10^(gex/20));
	gey=70;cgey=(10^(gey/20));
	ghx=30;cghx=(10^(ghx/20));
	ghy=30;cghy=(10^(ghy/20));
	ghz=30;cghz=(10^(ghz/20));

%-------------------------------------------------
%  Module for loading data
%-------------------------------------------------     
	load /home/dtrad/MT/datos.dat;
	tempdat= datos(1:ncut,1:5); % Cut the index greater than ncut
	datos=tempdat;
	clear tempdat;
	MMM=max(size(datos));

	if(inmax>MMM) 
	datos(MMM+1:inmax,:)=datos(MMM-1:-1:2*MMM-inmax,:);
	disp('The data will be wraped around');pause;
	end;
	
	ex=detrend(datos(inmin:inmax,1));
	ey=detrend(datos(inmin:inmax,2));
	hx=detrend(datos(inmin:inmax,3));
	hy=detrend(datos(inmin:inmax,4));
	hz=detrend(datos(inmin:inmax,5));

	clear datos;	
%--------------------------------------------------
%  Module for passband filtering (Fourier)
%--------------------------------------------------
	if(GAN=='y')
	ex=ex./cgex;
	ey=ey./cgey;
	hx=hx./cghx;
	hy=hy./cghy;
	hz=hz./cghz;
	end;


	if(FILT=='y')
	[bf,af]=butter(2,[0.01 0.9]);

	ex=filtfilt(bf,af,ex);
	ey=filtfilt(bf,af,ey);
	hx=filtfilt(bf,af,hx);
	hy=filtfilt(bf,af,hy);
	hz=filtfilt(bf,af,hz);

	end;
%--------------------------------------------------
%   Module for setting initial variables
%--------------------------------------------------
	n  = max(size(ex));
	D=log(n)/log(2);
	t= (1:n)./FS;	
 	QMF  = MakeONFilter(Family,Nfilt)
	if(GRAPHDATA=='y')
%--------------------------------------------------
% Plotting initial data.
%--------------------------------------------------
	subplot(221); 
	plot(t,ex)
	text=sprintf(' Observed ex(t) signal');
	title(text)
	
	subplot(222); 
	plot(t,ey)
	text=sprintf(' Observed ey(t) signal');
	title(text)

	subplot(223); 
	plot(t,hx)
	text=sprintf(' Observed hx(t) signal');
	title(text)

	subplot(224);
	plot(t,hy)
	text=sprintf(' Observed hy(t) signal');
	title(text)
	figure(gcf);pause,
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
	[exh,hxwh]=WaveShr1(ex,tipo,LS,QMF);
	[eyh,hywh]=WaveShr1(ey,tipo,LS,QMF);
	[hxh,hxwh]=WaveShr1(hx,tipo,LS,QMF);
	[hyh,hywh]=WaveShr1(hy,tipo,LS,QMF);
	[hzh,hzwh]=WaveShr1(hz,tipo,LS,QMF);

	elseif (filtro=='wpdeno')
	filtro
	[exh,bb,st] = wpdenoi1(ex,D,QMF);
	[eyh,bb,st] = wpdenoi1(ey,D,QMF);
	[hxh,bb,st] = wpdenoi1(hx,D,QMF);
	[hyh,bb,st] = wpdenoi1(hy,D,QMF);
	[hzh,bb,st] = wpdenoi1(hz,D,QMF);

	elseif (filtro=='mtfilt')
	filtro
    	if(FREQ=='y')
	
	windhan=hanning(164);
	wind=ones(size(ex));
	wind(1:82)=windhan(1:82);
	wind((8192-82+1):8192)=windhan(83:164);

	[mm nn]=size(wind);if(mm<nn) wind=wind';end;	

%	plot(wind),figure(gcf)

	 for i=1:1%Nwindows,
	  a=512*(i-1)+1;	
	  b=a+512-1;	 
	  a=1;b=8192; 

	  ex(a:b) = detrend(ex(a:b)).*wind;
	  ey(a:b) = detrend(ey(a:b)).*wind;
	  hx(a:b) = detrend(hx(a:b)).*wind;
	  hy(a:b) = detrend(hy(a:b)).*wind;
	  hz(a:b) = detrend(hz(a:b)).*wind;

	  	
	  exh(a:b) = fmtfiltro(ex(a:b),LS,QMF,option,tres,GRAPHS,GRAPHW,'ex','n');
	  eyh(a:b) = fmtfiltro(ey(a:b),LS,QMF,option,tres,GRAPHS,GRAPHW,'ey','n');
	  hxh(a:b) = fmtfiltro(hx(a:b),LS,QMF,option,tres,GRAPHS,GRAPHW,'hx','n');
	  hyh(a:b) = fmtfiltro(hy(a:b),LS,QMF,option,tres,GRAPHS,GRAPHW,'hy','n');
	  hzh(a:b) = fmtfiltro(hz(a:b),LS,QMF,option,tres,GRAPHS,GRAPHW,'hz','n');
	 end,
	else
	
	exh = mtfiltro(ex,LS,QMF,option,tres,GRAPHS,GRAPHW,'ex');
	eyh = mtfiltro(ey,LS,QMF,option,tres,GRAPHS,GRAPHW,'ey');
	hxh = mtfiltro(hx,LS,QMF,option,tres,GRAPHS,GRAPHW,'hx');
	hyh = mtfiltro(hy,LS,QMF,option,tres,GRAPHS,GRAPHW,'hy');
	hzh = mtfiltro(hz,LS,QMF,option,tres,GRAPHS,GRAPHW,'hz');
	end

	elseif (filtro=='multi1')
	filtro

    	exh = multimt1(ex,LS,QMF);
	eyh = multimt1(ey,LS,QMF);
	hxh = multimt1(hx,LS,QMF);
	hyh = multimt1(hy,LS,QMF);
	hzh = multimt1(hz,LS,QMF);

	end
%--------------------------------------------------------
	
	subplot(221); 
	plot(t,exh)
	text=sprintf(' Wavelet Reconstruction from ex filtered signal');
	title(text)
	
	subplot(222); 
	plot(t,eyh)
	text=sprintf(' Wavelet Reconstruction from ey filtered signal');
	title(text)

	subplot(223); 
	plot(t,hxh)
	text=sprintf(' Wavelet Reconstruction from hx filtered signal');
	title(text)

	subplot(224);
	plot(t,hyh)
	text=sprintf(' Wavelet Reconstruction from hy filtered signal');
	title(text)
	figure(gcf);pause

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
%   Module for auxiliar shaping
%----------------------------------------------------------------------
	[mmm,nnn]=size(ex);if(mmm<nnn) ex=ex';end;
	[mmm,nnn]=size(ey);if(mmm<nnn) ey=ey';end;
	[mmm,nnn]=size(hx);if(mmm<nnn) hx=hx';end;
	[mmm,nnn]=size(hy);if(mmm<nnn) hy=hy';end;
	[mmm,nnn]=size(hz);if(mmm<nnn) hz=hz';end;

	[mmm,nnn]=size(exh);if(mmm<nnn) exh=exh';end;
	[mmm,nnn]=size(eyh);if(mmm<nnn) eyh=eyh';end;
	[mmm,nnn]=size(hxh);if(mmm<nnn) hxh=hxh';end;
	[mmm,nnn]=size(hyh);if(mmm<nnn) hyh=hyh';end;
	[mmm,nnn]=size(hzh);if(mmm<nnn) hzh=hzh';end;
%-----------------------------------------------------------------------
% module for substracting the noise
%-----------------------------------------------------------------------
	pexh=ex-exh;
	peyh=ey-eyh;
	phxh=hx-hxh;
	phyh=hy-hyh;
	phzh=hz-hzh;
%------------------------------------------------------------------------
	subplot(221); 
	plot(t,pexh)
	text=sprintf('Removed noise from ex(t) signal');
	title(text)
	
	subplot(222); 
	plot(t,peyh)
	text=sprintf('Removed noise from ey(t) signal');
	title(text)

	subplot(223); 
	plot(t,phxh)
	text=sprintf('Removed noise from hx(t) signal');
	title(text)

	subplot(224);
	plot(t,phyh)
	text=sprintf('Removed noise from hy(t) signal');
	title(text)
	figure(gcf);pause

%	noise removed sfd	

	sfd=[pexh peyh phxh phyh phzh]; % removed noise

	clear pexh;
	clear peyh;
	clear phxh;
	clear phyh;
	clear phzh;

	Ndatos=min(MMM,inmax);
	if(GRAPHWIND=='y')		
	for i=1:Nwindows,
	a=512*(i-1)+1;
	b=a+512-1;
	ttex=sprintf('Noise removed:ex(t). Window %d',i);
	ttey=sprintf('Noise removed:ey(t). Window %d',i);
	tthx=sprintf('Noise removed:hx(t). Window %d',i);
	tthy=sprintf('Noise removed:hy(t). Window %d',i);
	MM=max(max(sfd(a:b,1:4)));NN=min(min(sfd(a:b,1:4)));
	subplot(221),plot(t(a:b),sfd(a:b,1)),title(ttex);xlabel('time(sec)')
	axis([t(a) t(b) NN MM]);
	subplot(222),plot(t(a:b),sfd(a:b,2)),title(ttey);xlabel('time(sec)')
	axis([t(a) t(b) NN MM]);
	subplot(223),plot(t(a:b),sfd(a:b,3)),title(tthx);xlabel('time(sec)')
	axis([t(a) t(b) NN MM]);
	subplot(224),plot(t(a:b),sfd(a:b,4)),title(tthy);xlabel('time(sec)')
	axis([t(a) t(b) NN MM]);
	figure(gcf),pause
	end,
	end,
%       clean signal : sf
	if (GAN=='y')
	exh=exh.*cgex;
	eyh=eyh.*cgey;
	hxh=hxh.*cghx;
	hyh=hyh.*cghy;
	hzh=hzh.*cghz;
	end;
	sf=[exh(1:Ndatos) eyh(1:Ndatos) hxh(1:Ndatos) hyh(1:Ndatos) hzh(1:Ndatos)];
	sf=[ex(1:Ndatos) ey(1:Ndatos) hx(1:Ndatos) hy(1:Ndatos) hz(1:Ndatos)];
%-----------------------------------------------------------------------
	
% Built October 04, 1996 
% Daniel Trad, CONICET.CRICYT. Mendoza. Argentina.

   
    
   

