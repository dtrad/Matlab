function [x]=deconv_trace(x,w)
% [x]=deconv_trace(x,w)
% Deconvolve seismic section x with wavelet w
% Daniel Trad
x=seis_shape(x);
[NT,NP]=size(x);
for ii=1:NP,
   y=deconv(x(:,ii),w);
   y=padzeros(y,NT);
   y=y(:);
   x(:,ii)=y;
end;
