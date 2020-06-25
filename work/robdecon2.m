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
nspikes=50;
optionw='rick'
%lfilter=lx;
lfilter=20;
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

cc=reflectivity2(c_rev,lx,nspikes);

cc(1:it(2)-1)=0;
[x,a]=rf2waves(cc,lx);
[cc_obt]=filter2coeff(a,lx-1);cc_obt(1)=0;
a_lev=levinson(x,lx-1);cc(1)=0;  % PEO from Matlab

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure,
subplot(221);linesad(tt,x_rev);title('(a) waves, zero offset');
subplot(222);linesad(tt,a_rev);title('(b) filter');
subplot(223);linesad(tt,c_rev);title('(c) reflection coeff (obtained)');
unif_scale(gcf,221,2,3);

figure,
subplot(221);linesad(tt,x);title('(a) waves, zero offset');
subplot(222);linesad(tt,a);title('(b) filter');
subplot(223);linesad(tt,cc);title('(c) reflection coeff (given)');
subplot(224);linesad(tt,cc_obt);title('(d) reflection coeff (obtained)');
unif_scale(gcf,221);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolution with wavelet 

fb=70;
lw=64;
wav=seismwave(optionw,lw,dt,fb);lw=length(wav);
wavmp=minphase(wav);
wav=wav(1:18);
figure,
subplot(221);plot(wav);title('(a) True wavelet')
subplot(222);plot(real(wavmp(1:18)));title('(b) MP wavelet')
subplot(223);isitmp2(wav);title('(c) True wavelet in the Z- domain')
subplot(224);isitmp2(wavmp(1:18));title('(d) MP wavelet in the Z-domain')

wav=normalize(wav);
%wav=real(minphase(wav));
%x=cc;
x(1)=1;
x_temp=convlim(wav,x,lx);

% Trace
x_conv=x_temp(1:lx);  
%isitmp(x_conv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part I: Spike deconvolution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predictive decon on convolved trace

[x_spike_deconv,peo_spike,wm,x_recovered]=spikedecon(x_conv,bigdiag,1,lfilter,'n');

figure,
subplot(221),linesad(tt,real(x));
title('(a)');xlabel('time (s)');ylabel('amplitude')
subplot(222),linesad(tt(1:lfilter),peo_spike(1:lfilter));title('(b)');
title('(b)');xlabel('lag (s)');ylabel('amplitude')
subplot(223),linesad(tt(1:lx),x_spike_deconv(1:lx));title('(c)')
title('(c)');xlabel('time (s)');ylabel('amplitude')
subplot(224),linesad(tt(1:lfilter),real(wm(1:lfilter)));title('(d)');
title('(d)');xlabel('time (s)');ylabel('amplitude')


if optionr=='y'
	% Recover the original coefficients
	[cc_obt_conv]=filter2coeff(peo_spike,lx-1);cc_obt_conv(1)=0;
	cc_obt_conv=normalize(cc_obt_conv).*cc(it(1));
	figure;
	subplot(221);linesad(real(wav));title('(a) Wavelet');
	subplot(222);linesad(tt,x);title('(a) Full Trace');
	subplot(223);linesad(x_conv);title('(b) Band limited trace');
	subplot(224);linesad(tt,peo_spike);title('(c) filter from band limited trace');
	%unif_scale(gcf,222);
	figure;	
	subplot(221);linesad(tt,cc);title('(a) reflection coeff (given)');
	subplot(222);linesad(tt,cc_obt);title('(b) reflection coeff from Claerbout');
	subplot(223);linesad(tt,cc_obt_conv);title('(c) reflection coeff from band lim trace');
	res=resta(cc,cc_obt_conv);
	subplot(224);linesad(tt,res);title('(d) Residuals r.c. :full trace - band lim. trace');
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
subplot(221),linesad(peo_gap);title(strcat('(a) peo-gap, delay=',int2str(delay)))
subplot(222),linesad(real(w_gap));title('(b) w_gap=peo-gap^-^1');
subplot(223),linesad(x_deconv_gap);title('(c) x deconvolved with gap deconvolution filter')
subplot(224),linesad(x_recovered_gap);title('(d) x recovered with g0')

x_gap_dec=convlim(w_gap(1:delay),x_spike_deconv);

%figure,
%subplot(221),linesad(x_gap_dec);title(strcat('x-gap-dec, delay=',int2str(delay)))
%subplot(222),linesad(wav(1:lw));title('wav')
%subplot(224),linesad(x_deconv_gap(1:delay+lw));title('x deconvolved with gap deconvolution filter')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase II: Model driven deconvolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gap deconvolution filter can be obtained as
% g=spike_filter*hm 
% where hm=head=wm(1:delay);
ref=zeros(1,2*lfilter);ref(1)=1;

for delay=5:5:25
   
   hm=wm(1:delay);  % Head: leading part of wm (wavelet without reverberation
   hm=normalize(hm);
   gap_filter=convlim(peo_spike,hm,lx);  %----> g=f*hm
   w_gap=deconv_freq_domain(ref,gap_filter);         %----> w=f^-1
%   lfilter=delay+1;
   x_deconv_gap2=convlim(x_conv,gap_filter(1:lfilter),lx);   %----> z=g*x
   %x_deconv_gap2=convlim(x_spike_deconv,hm,lx);				 %----> z=hm*y	

   figure,
	subplot(221),linesad(tt,gap_filter);title(strcat('peo-gap, delay=',int2str(delay)))
	subplot(222),linesad(tt(1:lx),w_gap(1:lx));title(strcat('W_M, delay=',int2str(delay)))
	subplot(223),linesad(tt,x_deconv_gap2);title('Deconvolved x(t)')
   subplot(224);plot(tt(1:length(hm)),hm);title('(a) h_M');
   %unif_scale(gcf,223,3,3);
	   
   %%%%%%%%%%%%%%%%%%%%%%%%
   % Remove the head hm
   x_rec=deconv_freq_domain(x_conv,hm);lx_r=length(x_rec);
   x_dec_gap_rec=deconv_freq_domain(x_deconv_gap2,hm);lx_d=length(x_dec_gap_rec);
   
%   figure,
%   subplot(221);plot(hm);title('(a) head:Minimum Phase counterpart of the wavelet');
%   subplot(222);plot(wav);title('(b) original wavelet');
%   subplot(223),linesad(x_rec(2:lx_r));title('(c) x deconvolved with the head')
% 	 subplot(224),linesad(x_dec_gap_rec(2:lx_d));title('(d) x deconv. with gap dec.filt.and head')
%   unif_scale(gcf,223,3,4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part III: Remove the true wavelet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delay=20;
   hm=wm(1:delay);  % Head: leading part of wm (wavelet without reverberation
   hm=normalize(hm);
   gap_filter=convlim(peo_spike,hm,lx);  %----> g=f*hm
   w_gap=deconv_freq_domain(ref,gap_filter);         %----> w=f^-1
   lfilter=delay+1;
   x_deconv_gap2=convlim(x_conv,gap_filter(1:lfilter),lx);   %----> z=g*x
   %x_deconv_gap2=convlim(x_spike_deconv,hm,lx);				 %----> z=hm*y	

   figure,
   subplot(221),linesad(tt,gap_filter);title(strcat('gap filter, delay=',int2str(delay)))
   subplot(222),linesad(tt(1:lx),w_gap(1:lx));title(strcat('W_M, delay=',int2str(delay)))
	subplot(223),linesad(tt,x_deconv_gap2);title('x(t) deconvolved')
   subplot(224);plot(tt(1:length(hm)),hm);title('(a) h_M');
   %unif_scale(gcf,223,3,4);
   
   %%%%%%%%%%%%%%%%%%%%%%%%
   % Remove the head hm
   x_rec=deconv_freq_domain(x_conv,hm);lx_r=length(x_rec);
   x_dec_gap_rec=deconv_freq_domain(x_deconv_gap2,hm);lx_d=length(x_dec_gap_rec);
   
   %figure,
   %subplot(221);plot(hm);title('(a) head:Minimum Phase counterpart of the wavelet');
   %subplot(222);plot(wav);title('(b) original wavelet');
	%subplot(223),linesad(x_rec(2:lx_r));title('(c) x deconvolved with the head')
	%subplot(224),linesad(x_dec_gap_rec(2:lx_d));title('(d) x deconv. with gap dec.filt.and head')
   %unif_scale(gcf,223,3,4);
   wav0=padzeros(wav,lx);
   p=deconv_freq_domain(wav0,hm);lp=length(p);
   %pinv=deconv(ref,p);
   pinv=p(lp:-1:1);
   figure,
   subplot(221);plot(p);title('(a) phase shift filter p');
   subplot(222);plot(pinv);title('(b) inverse phase shift filter p')
   
   wavr=convlim(pinv,wav,lx);
   wavrr=shift_time(wavr,0.2);
   %subplot(223);plot(wavrr);title('(c) p^-^1*wav=MP wav'); 
  
  % Remove the true wavelet
   x_spike_deconv_mp=convlim(pinv,x_spike_deconv); 	% Remove the phase shift filter
   x_spike_deconv=deconv_freq_domain(x_spike_deconv,hm);   			% Remove the MP wavelet
   x_spike_deconv_mp=shift_time(x_spike_deconv_mp); 	% Shift time origin
   
   mp_dec_gap_rec_mp=convlim(pinv,x_dec_gap_rec); 		% Remove the phase shift filter
   mp_dec_gap_rec_mp=shift_time(mp_dec_gap_rec_mp);	% Shift time origin
   
   lx_r=length(x_spike_deconv_mp);
   lx_t=length(mp_dec_gap_rec_mp);
   
   temp1=max(abs(x(80:lx)));
   temp2=max(abs(mp_dec_gap_rec_mp(100:lx_t)));
   temp3=temp1./temp2;
   mp_dec_gap_rec_mp=mp_dec_gap_rec_mp*temp3;
      
   figure,
	subplot(221),linesad(x_spike_deconv_mp(1:lx_r));title('(a) spike deconv. trace MP ')
	subplot(222),linesad(mp_dec_gap_rec_mp(2:lx_t));title('(b) gap deconv. trace MP')
   subplot(223);linesad(real(x(2:lx)));title('(c) original trace');
      
   unif_scale(gcf,223,2,3);
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Predictive decon on convolved trace MP
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Remove the true wavelet
   
   x_conv_mp=convlim(pinv,x_conv);
   x_rec_mp=deconv_freq_domain(x_conv_mp,hm);
   x_rec_mp=shift_time(x_rec_mp);
   x_rec_mp=padzeros(x_rec_mp,lx);
   
   %lx_r=length(x_rec_mp);
   %figure,
	%subplot(223),linesad(x_conv_mp(2:lx_r));title('MP trace')
   
	[x_spike_deconv_mp,peo_spike_mp,wm_mp,x_recovered_mp]=spikedecon(x_rec_mp,bigdiag,1,lx,'y');

	figure,
	subplot(221),linesad(real(wm_mp));title('(a) wm=f^-^1');
	subplot(222),linesad(peo_spike_mp);title('(b) spike deconvolution filter');
	subplot(223),linesad(x_spike_deconv_mp);title('(c) x dec. with sp. dec. filter')
	subplot(224),linesad(x_recovered_mp);title('(d) Convolution of de-spike x and wm')
   
   
   

   
   