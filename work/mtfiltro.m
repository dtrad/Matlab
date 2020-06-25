% 	Function for filtering MT data with filtror.
%
% 	Called by mtwav03.m
%       Call to filtror.m
%
%	Input Arguments
%        dirty  = noisy signal
%        L      = coarser resolution
%        QMF    = wavelet filter
%        option = 'Huber', 'Thomson'
%        tres   = Number of standar deviations to consider acceptable.
%
%	Output Argument
%        clean  = denoised signal 
%
%       Daniel Trad- 10-10-96

	function [clean]=mtfiltro(dirty,L,QMF,option,tres,GRAPHS,GRAPHW,chan)
	
	wdirty = FWT_PO(dirty,L,QMF);
        wclean=filtror2(wdirty,L,tres,option,GRAPHS,chan);	
	clean =IWT_PO(wclean,L,QMF);

	if(GRAPHW=='y')
	subplot(111); plotwavecoeff(wclean-wdirty,L,0);
	xlabel('time');ylabel('scale')
	text=sprintf('WT: Removed %s noise',chan);title(text);
	figure(gcf);%pause
	subplot(111); plotmra(wclean-wdirty,L,0,QMF);
	xlabel('time');ylabel('scale')
	text=sprintf('MRA- Removed %s noise',chan);title(text);
	figure(gcf);%pause

	subplot(111); plotwavecoeff(wclean,L,0);
	xlabel('time');ylabel('scale')
	text=sprintf('WT: Denoised %s signal',chan);title(text);
	figure(gcf);%pause
	subplot(111); plotmra(wclean,L,0,QMF);
	xlabel('time');ylabel('scale')
	text=sprintf('MRA- Denoised %s signal',chan);title(text);
	figure(gcf);%pause
	end
