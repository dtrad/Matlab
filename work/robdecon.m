% Program Robdecon
% Zero offset deconvolution
% Robinson's example: Non minimum phase wavelet with 1st order reverberation 
% Theory:Robinson: Model Driven deconvolution
% Geophysics Vol63-No2 pag.1199-1206
% GEOP 520: Direct Studies
% Daniel Trad: UBC. 
echo on
clear;
myclean;
dt=0.004;
lx=256;
option=1;
bigdiag=1.0;
nspikes=150;
optionw='rick'
lfilter=lx;
optionr='n';
optiongap1='n';
%optionw='minr'
%optionw='kais'
%optionw='no_w';

% First order reverberation filter
t=[0 .08];  
rf=[1 .5];
n=length(rf);

it=round(t./dt+1);
a=zeros(1,lx); 
a(it)=rf;
tt=0:dt:((lx-1)*dt);

[c_rev]=filter2coeff(a,lx-1);
[x_rev,a_rev]=rf2waves(c_rev,lx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full reflectivity series and waves

cc=reflectivity2(x_rev,lx,nspikes);
cc(1:it(1)-1)=0;
[x,a]=rf2waves(cc,lx);
[cc_obt]=filter2coeff(a,lx-1);cc_obt(1)=0;
a_lev=levinson(x,lx-1);cc(1)=0;  % PEO from Matlab

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure,
subplot(221);linesad(tt,x_rev);title('waves, zero offset');
subplot(222);linesad(tt,a_rev);title('filter');
subplot(224);linesad(tt,c_rev);title('reflection coeff (obtained)');
unif_scale(gcf,221);

figure,
subplot(221);linesad(tt,x);title('waves, zero offset');
subplot(222);linesad(tt,a);title('filter');
subplot(223);linesad(tt,cc);title('reflection coeff (given)');
subplot(224);linesad(tt,cc_obt);title('reflection coeff (obtained)');
unif_scale(gcf,221);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolution with wavelet 

fb=70;
lw=64;
wav=seismwave(optionw,lw,dt,fb);lw=length(wav);
wav=wav(1:18);
figure,isitmp2(wav);
wav=normalize(wav);
%wav=real(minphase(wav));
x=cc;
x(1)=1;
x_temp=convlim(wav,x,lx);

% Trace
x_conv=x_temp(1:lx);  
%isitmp(x_conv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predictive decon on convolved trace

[x_spike_deconv,peo_spike,wm,x_recovered]=spikedecon(x_conv,bigdiag,1,lfilter,'y');

figure,
subplot(221),linesad(real(wm));title('wm=f^-^1');
subplot(222),linesad(peo_spike);title('spike deconvolution filter');
subplot(223),linesad(x_spike_deconv);title('x dec. with sp. dec. filter')
subplot(224),linesad(x_recovered);title('Convolution of de-spike x and wm')

if optionr=='y'
	% Recover the original coefficients
	[cc_obt_conv]=filter2coeff(peo_spike,lx-1);cc_obt_conv(1)=0;
	cc_obt_conv=normalize(cc_obt_conv).*cc(it(1));
	figure;
	subplot(221);linesad(wav);
	subplot(222);linesad(tt,x);title('Full Trace');
	subplot(223);linesad(x_conv);title('Band limited trace');
	subplot(224);linesad(tt,peo_spike);title('filter from band limited trace');
	%unif_scale(gcf,222);
	figure;	
	subplot(221);linesad(tt,cc);title('reflection coeff (given)');
	subplot(222);linesad(tt,cc_obt);title('reflection coeff from Claerbout');
	subplot(223);linesad(tt,cc_obt_conv);title('reflection coeff from band lim trace');
	res=resta(cc,cc_obt_conv);
	subplot(224);linesad(tt,res);title('Residuals r.c. :full trace - band lim. trace');
	%unif_scale(gcf,221);
end;
if optiongap1=='y'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different gaps
%bigdiag=1;
figure
for ii=1:3;
   delay=2*ii+13;
   [x_deconv_gap,peo_gap,w_gap,x_recovered_gap]=spikedecon(x_conv,bigdiag,delay,lfilter);
   subplot(229+2*ii);linesad(peo_gap);title(strcat('peo-gap, delay=',int2str(delay)))
	subplot(230+2*ii),linesad(real(w_gap));title('w=peo-gap^-^1');
end;

delay=20;
[x_deconv_gap,peo_gap,w_gap,x_recovered_gap]=spikedecon(x_conv,bigdiag,delay,lfilter);

figure,
subplot(221),linesad(peo_gap);title(strcat('peo-gap, delay=',int2str(delay)))
subplot(222),linesad(real(w_gap));title('w_gap=peo-gap^-^1');
subplot(223),linesad(x_deconv_gap);title('x deconvolved with gap deconvolution filter')
subplot(224),linesad(x_recovered_gap);title('x recovered with g0')

x_gap_dec=convlim(w_gap(1:delay),x_spike_deconv);

%figure,
%subplot(221),linesad(x_gap_dec);title(strcat('x-gap-dec, delay=',int2str(delay)))
%subplot(222),linesad(wav(1:lw));title('wav')
%subplot(224),linesad(x_deconv_gap(1:delay+lw));title('x deconvolved with gap deconvolution filter')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gap deconvolution filter can be obtained as
% g=spike_filter*hm 
% where hm=head=wm(1:delay);
ref=zeros(1,2*lfilter);ref(1)=1;

for delay=18:22
   
   hm=wm(1:delay);  % Head: leading part of wm (wavelet without reverberation
   hm=normalize(hm);
   gap_filter=convlim(peo_spike,hm,lx);  %----> g=f*hm
   w_gap=deconv(ref,gap_filter);         %----> w=f^-1
%   lfilter=delay+1;
   x_deconv_gap2=convlim(x_conv,gap_filter(1:lfilter),lx);   %----> z=g*x
   %x_deconv_gap2=convlim(x_spike_deconv,hm,lx);				 %----> z=hm*y	

   figure,
	subplot(221),linesad(w_gap);title(strcat('peo-gap, delay=',int2str(delay)))
	subplot(222),linesad(gap_filter);title(strcat('gap_filter, delay=',int2str(delay)))
	subplot(223),linesad(x_conv);title('x')
	subplot(224),linesad(x_deconv_gap2);title('x deconvolved with gap deconvolution filter 2')
   unif_scale(gcf,223,3,4);
   
   %%%%%%%%%%%%%%%%%%%%%%%%
   % Remove the head hm
   x_rec=deconv(x_conv,hm);lx_r=length(x_rec);
   x_dec_gap_rec=deconv(x_deconv_gap2,hm);lx_d=length(x_dec_gap_rec);
   
   figure,
   subplot(221);plot(hm);title('head');
   subplot(222);plot(wav);title('original wavelet');
	subplot(223),linesad(x_rec(2:lx_r));title('x')
	subplot(224),linesad(x_dec_gap_rec(2:lx_d));title('x deconvolved with gap deconvolution filter 2')
	   
	unif_scale(gcf,223,3,4);
end