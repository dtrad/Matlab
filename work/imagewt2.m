ingrid = ReadImage('Daubechies');
%       
        figure(1)
	subplot(111)
	GrayImage(ingrid)
	title('Ingrid Daubechies');

	qmf = MakeONFilter('Coiflet',2);
	wingrid = FWT2_PO(ingrid,3,qmf);
%
	zmat = abs(wingrid);
	figure(2)
	AutoImage(zmat);colorbar
	title('Wavelet Transform of Ingrid Daubechies');
	
%  Investigate Sparsity in the Wavelet Transform of Daubechies.
%
	wcsort = sort(abs(wingrid(:)));
	wcerr  = cumsum(wcsort.^2);
	wcerr  = flipud(wcerr);
%
%  Sparsify Image
%
	wthresh = wcsort(floor(.95*65536));
	cw_ingrid = wingrid .* (abs(wingrid) > wthresh);
	[i,j,s] = find(cw_ingrid);
	sp_ingrid = sparse(i,j,s,256,256);
	figure(3)
	spy(sp_ingrid)
	title('Nonzero Pattern in Sparsification of WT[Daubechies]')


%  Reconstruct Daubechies from 5% of her coefficients.
%
	icw_ingrid = IWT2_PO(cw_ingrid,3,qmf);
	figure(4)
	AutoImage(icw_ingrid);
	title('95% Wavelet Co/Dec of Daubechies');
    
    

	