function [y,res]=wav1dfilter(x,tres,LS)
% wavelet filter 
% Program for filtering time series contaminated by random noise
% example:
% [y,res]=wav1dfilter(x+n,1,nextpow2(length(x+n))-2);

% calls: wavfilt.m---> subcalls: wavfilter.m
% need : wavelab
%
% Daniel Trad - 26-10-96 - Updated June 17, 2002
% Default parameters are reasonable but they can be changed
% in the first part of the code	


	if (nargin<3) LS=nextpow2(length(x+n))-2;end %coarser level
	if (nargin<2) tres=1.0;end %treshold

%------------------------------------------------
%  change parameters here if you need.
%------------------------------------------------  
	option=1 %Huber (mtfilt)
% 	option=2 %Thomson (mtfilt)
%	filtro='Shrink'
%	filtro='wpdeno'
	filtro='mtfilt'
%	filtro='multi1'
%	Family='Haar',Nfilt=2;
% 	Family='Daubechies', Nfilt=4;
	Family='Symmlet', Nfilt=8;
% 	Family='Coiflet', Nfilt=3;

	FS=0.004 %sample frequency (only for plots)
	
	tipo='MAD';
	% The series need to have a power of two
	% Can be greater than size of data.
	GRAPHS='y';  % Graphics of Wavelet coefficient series
	GRAPHW='n';  % Graphics of Wavelet Coefficients
	fign=1; % first figure number

%-----------------------------------------------------------------
	lx=length(x);
	lx2=2^(nextpow2(x)); % Need length of data equal to power of two
	% Pad with the data themselves to reach power of two
	if(lx2>lx) 
	  x(lx+1:lx2)=x(lx-1:-1:2*lx-lx2,:);
	  disp('The data will be wraped around');pause;
	end;

	t=1:length(x);t=t*FS;
	figure(fign);fign=fign+1;
	plot(t,x);
	text=sprintf(' Observed x(t) signal');
	title(text)
	
	x=detrend(x);

        n  = max(size(x));
	D=log(n)/log(2);
	t= (1:n)./FS;	
 	QMF  = MakeONFilter(Family,Nfilt)

	if(filtro=='Shrink')
	filtro
	[y,yw]=WaveShr1(x,tipo,LS,QMF);

	elseif (filtro=='wpdeno')
	filtro
	[y,bb,st] = wpdenoi1(x,D,QMF);

	elseif (filtro=='mtfilt')
	filtro
	y = wavfilt(x,LS,QMF,option,tres,GRAPHS,GRAPHW,'x');

	elseif (filtro=='multi1')
	filtro
    	y = multimt1(x,LS,QMF);
	end
%--------------------------------------------------------
        figure(fign);fign=fign+1;
        subplot(311);
	plot(t,x)
	text=sprintf('Original signal ');
	title(text);
	V=axis;

	
        subplot(312);
	plot(t,y)
	text=sprintf(' Wavelet Reconstruction ');
	title(text)
	axis(V);
	
	res=x-y;
	
        subplot(313);
	plot(t,res)
	text=sprintf('Removed noise ');
	title(text)
	axis(V);
	
%-----------------------------------------------------------------------
	
% Built October 04, 1996 
% Daniel Trad, CONICET.CRICYT. Mendoza. Argentina.
% Adapted for single time series with random noise June 17, 2002.
% Daniel Trad, UBC. 

   
    
   


