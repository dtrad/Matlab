function [wmest,wcmest,twcmest]=applywt2(mest,nt,nh,perc)
% generate the wavelet
L = 1;

qmf = MakeONFilter('Battle',3);

if (nt==nh)
  display('Fast 2d FWT');
  wcmest=FWT2_PO(mest,L,qmf);
else
  % do the 2-D wavelet transform in two passes
  for itr=1:nh
    wcmest(:,itr) =  FWT_PO(mest(:,itr),L,qmf);
  end

  for is=1:nt
    wcmest(is,:) =  FWT_PO(wcmest(is,:),L,qmf);
  end

  % compute the noise level
end

% do the thresholding 
%twcmest =  HardThresh(wcmest,T);

[T]=threshold(wcmest,perc);
twcmest=wcmest;
twcmest(abs(wcmest)<T)=0;
%twcmest=wcmest;
% do the reconstruction

if (nt==nh)
  display('Fast 2d IWT');
  wmest=IWT2_PO(twcmest,L,qmf);
else
  for itr=1:nh,
    wmest(:,itr) =  IWT_PO(twcmest(:,itr),L,qmf);
  end

  for is=1:nt
    wmest(is,:) =  IWT_PO(wmest(is,:),L,qmf);
  end
end



