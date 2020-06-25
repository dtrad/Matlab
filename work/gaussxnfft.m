function [gaussx] = nfft_lookup(lenint, nlookup, iover)


alpha = (2.0 - 1.0/iover)/lenint
beta = pi*alpha

ipos= 0:nloop-1;
xpos = ipos/(nlookup-1);
is = -lenint/2+1:lenint/2;

iss = is(:)*ones(1,length(xpos));
iposs = ones(length(is),1)*xpos(:).';

gaussx = exp(-beta*(xposs - iss)^2);

return;
