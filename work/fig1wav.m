% toon0131 -- Scale Families of Wavelets
%
%  Show the Symmlet 8 wavelet at various scales and locations
%
	clf; subplot(111);
	posarray = [ 3 2 ; 3  5; 4 8; 5 13; 6 21; 6 32 ; 6 43; 7 95 ];
	sz = size(posarray);
	nr = sz(1);
	n  = 1024;
	w = zeros(1,n);
	t = (.5:(n-.5)) ./n;
%
	LockAxes([0 1 0 (nr+1)]); 
	%title('')
        xlabel('time')
        ylabel('scale');

	for iter = 1:nr,
		j = posarray(iter,1);
		k = posarray(iter,2);
		w = MakeWavelet(j,k,'Daubechies',4,'Mother',1024);
		plot(t,(iter)+ 3*w);
		txt = sprintf('(%1.0f,%2.0f)',j,k);
		text(.87,(iter)+.275,txt);
	end

	UnlockAxes;
    
    
%   
% Part of WaveLab Version .701
% Built Tuesday, January 30, 1996 8:25:59 PM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@playfair.stanford.edu
%   
    
