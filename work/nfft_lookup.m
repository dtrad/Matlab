function [gaussx] = nfft_lookup(lenint, nlookup, iover)


alpha = (2.0 - 1.0/iover)/lenint
beta = pi*alpha

ipos= 0:nlookup-1;
xpos = ipos/(nlookup-1);
is = -lenint/2+1:lenint/2;

iss = is(:)*ones(1,length(xpos));
xposs = ones(length(is),1)*xpos(:).';

%iss(1:10,1:10)

%xposs(1:10,1:10)

gaussx = exp(-beta*(xposs - iss).^2);

return;
