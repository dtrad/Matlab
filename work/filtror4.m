function [bex] = filtror4(wex,LS,treshold,option,GRAPHS,PART)
% 	Function for filtering frequency MT data in wavelet domain.
%	It is computed a treshold and scale factor for each level (i) and 
%       a weight for each coeficient along each level.
%	
% 	Called by mtfiltro.m
%       
%
%	Input Arguments
%        wex  = real or imaginary part of noisy signal in wavelet domain
%        LS      = coarser resolution
%        treshold   = Number of standar deviations to consider acceptable.
%        option = 'Huber', 'Thomson'
%	 GRAPHS = option for plotting
%	
%	Output Argument
%        bex  = real or imag. denoised signal in wavelet domain
%       Daniel Trad- 10-10-96

	bex=wex;

	w=ones(size(wex));
	if (option==1) turn2=treshold;end
	[n,J]=dyadlength(wex);
%-----------------------------------------------
% Module for adaptative scale filtering
%-----------------------------------------------	
% The filter computes a weight scale by scale and apply on
% the real and imaginary parts of the MT signal avoiding
% changes in phase.
% i defines the level
% j defines the traslation

	for i=J-1:-1:LS,
		clear a;
		a=(wex(dyad(i)));
		tolx=40;
		if (option==2|option==4) turn2=(J/J)*(2*log(2*max(size(a))))^0.5;end
		scale=1/0.44845*median(abs(a-median(a)));
		for j=min(dyad(i)):max(dyad(i)),
			srx=abs(wex(j))./scale;
			if (option==1)	%Huber	
			   if (srx < treshold)
				w(j)=1;
		    	   elseif (srx > treshold) & (srx < (3*treshold))
				w(j)=(treshold./srx);
			   else
			        w(j)=0; (treshold./srx)^2;	
			   end
			elseif (option==2|option==4)  % Thomson
			   expo=turn2*(srx-turn2);		
			   if(expo > tolx) 
				w(j)=0;
			   else    
			 	w(j)=exp(-exp(expo));
			   end
			end
		end
		if(option==4)
			clear a;
			dd=dyad(i);dyadcent=dd(2^(i-2)+1:2^(i-2)+2^(i-1));
			a=(wex(dyadcent));
			tolx=40;
			scalecent=1/0.44845*median(abs(a-median(a)));
 			turn2cent=0.5*(J/J)*(2*log(2*max(size(a))))^0.5;
			for j=min(dyadcent):max(dyadcent),
				srx=abs(wex(j))./scalecent;
				expo=turn2cent*(srx-turn2cent);		
				if(expo > tolx) 
					w(j)=0;
				else    
			 		w(j)=exp(-exp(expo));
				end
			end
		end
		if(GRAPHS=='y')
		  subplot(111),
		  t=1:max(size(dyad(i)));
		  turn2m=turn2*ones(size(dyad(i)));
		  plot(t,wex(dyad(i))./scale,t,turn2m,t,-1.*turn2m,t,wex(dyad(i)).*w(dyad(i))./scale),
		  text=sprintf('level %d',i);title(text);
		  xlabel('traslation');ylabel(PART),figure(gcf);pause
		end
	end
	
	%plot(w),title('Weights');figure(gcf);
	bex=wex.*(w);

