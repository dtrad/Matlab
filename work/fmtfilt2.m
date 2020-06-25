% 	Function for filtering freq MT data with filtror4 in DWT.
%	It differs from fmtfiltr.m in that apply the weights independently 
%       to Real and Imaginary parts. (It uses filtror4.m instead of filtror3.m)
%
% 	Called by mtwav04.m
%       Call to filtror4.m
%
%	Input Arguments
%        dirty  = noisy signal
%        L      = coarser resolution
%        QMF    = wavelet filter
%        option = 'Huber', 'Thomson'
%        tres   = Number of standar deviations to consider acceptable.
%	 GRAPHS = Option for plotting Series of Wavelet Coeficients	
%	 GRAPHW = Option for plotting Wavelet Coeficients (Bi-parametric)	
%	Output Argument
%        clean  = denoised signal 
%
%       Daniel Trad- 10-10-96

	function [clean]=fmtfilt2(dirty,L,QMF,option,tres,GRAPHS,GRAPHW,chan,PLOTSP)
	
	dirtyr=real(fft(dirty));
	dirtyi=imag(fft(dirty));
	dirty=(dirtyr.^2+dirtyi.^2).^0.5;
	
	wdirty = FWT_PO(dirty,L,QMF);
	wdirtyr = FWT_PO(dirtyr,L,QMF);
	wdirtyi = FWT_PO(dirtyi,L,QMF);

        [wcleanr]=filtror4(wdirtyr,L,tres,option,GRAPHS,chan);	
        [wcleani]=filtror4(wdirtyi,L,tres,option,GRAPHS,chan);	

	cleanr =IWT_PO(wcleanr,L,QMF);
	cleani =IWT_PO(wcleani,L,QMF);
	cleant=cleanr+cleani.*((-1)^0.5);

	cleanm=(cleanr.^2+cleani.^2).^0.5;
	clean=real(ifft(cleant));

	if (PLOTSP=='y')
	tam=max(size(cleanm));
	ff=(0:tam/2-1).*500./tam;

	subplot(211) 
	plot(ff,dirty(1:tam/2));xlabel('Frequency(Hz)');ylabel('Amplitude(mV)');
	subplot(212) 
	plot(ff,cleanm(1:tam/2));xlabel('Frequency(Hz)');ylabel('Amplitude(mV)');
	figure(gcf);pause;

	subplot(111)
	plot(ff,dirty(1:tam/2)-cleanm(1:tam/2));xlabel('Frequency(Hz)');ylabel('Amplitude(mV)');
	figure(gcf);pause;

	end,

	if(GRAPHW=='y')
	subplot(111); PlotWaveCoeff(wcleanr-wdirtyr,L,0);
	xlabel('time');ylabel('scale')
	text=sprintf('WT: Removed Real %s noise',chan);title(text);

	%subplot(111); PlotMRA_MT(wcleanr-wdirtyr,L,0,QMF);
	%xlabel('time');ylabel('scale')
	%text=sprintf('MRA- Removed Real %s noise',chan);title(text);
	%figure(gcf);pause

	subplot(111); PlotWaveCoeff(wcleanr,L,0);
	xlabel('time');ylabel('scale')
	text=sprintf('WT: Denoised Real %s signal',chan);title(text);
	figure(gcf);pause

	%subplot(111); PlotMRA_MT(wcleanr,L,0,QMF);
	%xlabel('time');ylabel('scale')
	%text=sprintf('MRA- Denoised Real %s signal',chan);title(text);
	%figure(gcf);pause
	end
