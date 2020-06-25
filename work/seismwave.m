function [wav,lw]=seismwave(optionw,lo,dt,fb)
% [w]=seismwave(optionw,lx,dt,fb)
% Select the seismic wavelet to use
% defaults: optionw='rick';
%				lx=50; (length)
%           dt=0.004;
%           fb=25;    (Central freq. for ricker wavelet) 
% Daniel Trad UBC. 30-07-98
if nargin<4 fb=25;end
if nargin<3 dt=0.004;end
if nargin<2 lo=50;end
if nargin<1 optionw='rick';end
optionw=optionw(1:4);
if optionw=='rick'
   wav=rickerm(fb,dt);
   wav=padzeros(wav,lo);
   wav=wav(1:lo);
   lw=length(wav);
	%figure,plot(wav)
elseif optionw=='minr'
   wavr=rickerm(fb,dt);
   wav=padzeros(wavr,lo);
   wav=minphase(wav);lw=length(wav);
   figure,
   subplot(221),plot(real(wav));title('MP wavelet');
   subplot(222),plot(imag(wav));title('Imag Part ~ 0');
   subplot(223),plot(wavr);title('Ricker NMP');

   wav=real(wav);
elseif optionw=='kais'
   load kaiserspike;lw=length(temp);
   inipoint=50;wav=temp(inipoint:lw-inipoint+1);
   figure,plot(wav)
elseif optionw=='spik'
   wav=1;lw=1;    
end

