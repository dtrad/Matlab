  	function [bex] = filtror2(wex,LS,treshold,option,GRAPHS,chan)
% 	Function for filtering MT data in wavelet domain.
%	It is computed a treshold and scale factor for each level (i) and 
%       a weight for each coeficient along each level.
%	
% 	Called by mtfiltro.m
%       
%
%	Input Arguments
%        wex  = noisy signal in wavelet domain
%        LS   = coarser resolution
%        treshold = Number of standar deviations to consider acceptable.
%        option = 'Huber', 'Thomson'
%
%	Output Argument
%        bex  = denoised signal in wavelet domain
%
%       Daniel Trad- 10-10-96
		
	bex=wex;
	w=ones(size(wex));

	[n,J]=dyadlength(wex);
%-----------------------------------------------
% Module for adaptative scale filtering
%-----------------------------------------------	
% The filter apply downweighting scale by scale.
% i defines the level
% j defines the traslation

	for i=J-1:-1:LS, 
		clear a;
		a=(wex(dyad(i)));
		tolx=40;
		if (option==1) turn2=(J/i) * treshold;end;
		if (option==2) turn2=(J/i)*(2*log(max(size(a))))^0.5;end
		
		scale=1.4826*median(abs(a-median(a)));

		for j=min(dyad(i)):max(dyad(i)), 
			srx=abs(wex(j))./scale;
			if (option==1)	%Huber	
			   if (srx < turn2)
				w(j)=1;
		    	   elseif (srx > turn2) & (srx < (3*turn2))
				w(j)=(turn2./srx);
			   else
			        w(j)=0; (turn2./srx)^2;	
			   end
			elseif (option==2)  % Thomson
			   expo=turn2*(srx-turn2);		
			   if(expo > tolx) 
				w(j)=0;
			   else    
			 	w(j)=exp(-exp(expo));
			   end
			end
		end

% 		Plotting each level
		if(GRAPHS=='y')
		subplot(111),
		t=1:max(size(dyad(i)));
		turn2m=turn2*ones(size(dyad(i)));
		plot(t,wex(dyad(i))./scale,t,turn2m,t,-1.*turn2m,t,wex(dyad(i)).*w(dyad(i))./scale-120),
		text=sprintf('%s(%d)',chan,i);title(text);
		xlabel('Translation');ylabel('Amplitude'),
		axis([0 2200 -130 80]);
		figure(gcf);pause;
		end
	end

	%subplot(111)
   %plot(w),title('Weights on DWT');figure(gcf);
   exo=bex;
	bex=wex.*(w);
   if chan=='ex' save wex.mat bex; end
   if chan=='ex' save exo.mat exo; end

