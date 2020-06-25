% MFWork01: Multifractal Workout: CWT Analysis of Cantor Signal
clear 
close all
	load medfilt
	% modifiable parameters
	[N N]=size(medfilt);;     % signal length; fairly large
	nvoice = 9;  % Following JS Bach, Well-Tempered Klavier
	t=1:N;

	% create a Brownian
%	CantorMeasure = MakeFractal(N,3,'Deterministic',[.5 0 .5]);
%	Devil  = cumsum(CantorMeasure); t = (.5:(N-.5))./N;
        Devil=medfilt(:,90);
figure(1); plot(t, Devil); title(sprintf(' Devil Staircase Signal'));
	
	% make CWT
	% Devil = Devil - Devil(N) .*t;	
	Devil_cwt = CWT(Devil,nvoice,'Sombrero');
	sz = size(Devil_cwt); nscale = sz(2);

	% display CWT
	figure(2); ImageCWT(Devil_cwt,'Individual','hot');
	title('CWT')
	
	% Build Maxima Map
	Devil_maxmap = WTMM(Devil_cwt);
	% display maxmap
	figure(3); ImageWTMM(Devil_maxmap)
	
	% Identify Ridges
	[skm,skp,skl] = BuildSkelMapFast(Devil_maxmap);

	% display Ridges
	figure(4); PlotSkelMap(N,nscale,skm,skp,skl);
	
	% RidgePlots
	ridgelist = [ 6, 7, 13, 29, 76, 115];
	figure(5); PlotRidges(ridgelist,Devil_cwt,skm,skp,skl);
	
	% PruneRidges
	[skellist,skelptr,skellen] = PruneSkelMap(Devil_cwt,.001,1,skm,skp,skl);
	figure(6); PlotSkelMap(N,nscale,skellist,skelptr,skellen);
    
 
	return;

%	threshold=max(max(
	  
	  [n(ii) scale(ii)]=find((Devil_cwt)==max(max(Devil_cwt(:,: ...
							    ))));
	  

	  [n m]=find(Devil_cwt==min(min((Devil_cwt(150:250,:))))); 
	
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    
