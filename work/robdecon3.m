% Program Robdecon
% Zero offset deconvolution: 
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
bigdiag=1;
nspikes=15;
optionw='rick'
%optionw='minr'
%optionw='kais'
%optionw='no_w';
% First order reverberation filter
t=[0.08 0.16];  
rf=[1 .25];rf=rf;
n=length(rf);

it=round(t./dt+1);
a=zeros(1,lx); 
a(it)=rf;
tt=0:dt:((lx-1)*dt);
[c_rev]=filter2coeff(a,lx-1);
[R_rev,a_rev]=rf2waves(c_rev,lx);

cc=reflectivity2(c_rev,lx,nspikes);
cc(1:it(1)-1)=0;
[R,a]=rf2waves(cc,lx);
[cc_obt]=filter2coeff(a,lx-1);cc_obt(1)=0;
a_lev=levinson(R,lx-1);cc(1)=0;  % PEO from Matlab

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure,
subplot(221);linesad(tt,R_rev);title('waves, zero offset');
subplot(222);linesad(tt,a_rev);title('filter');
subplot(224);linesad(tt,c_rev);title('reflection coeff (obtained)');
unif_scale(gcf,221);

figure,
subplot(221);linesad(tt,R);title('waves, zero offset');
subplot(222);linesad(tt,a);title('filter');
subplot(223);linesad(tt,cc);title('reflection coeff (given)');
subplot(224);linesad(tt,cc_obt);title('reflection coeff (obtained)');
unif_scale(gcf,221);


fb=25;
lw=30;
wav=seismwave('spike',lw,dt,fb);lw=length(wav);
wav=normalize(wav);
R(1)=1;

R_temp=convlim(wav,R,lx);

% Trace
R_conv=R_temp(1:lx);  
%isitmp(R_conv);
ref=zeros(1,2*lx);ref(1)=1;

[x_spike_deconv,peo_spike,wm,x_recovered]=spikedecon(R_conv,bigdiag,1,'y');

figure,
subplot(221),linesad(real(wm));title('wm=f^-^1');
subplot(222),linesad(peo_spike);title('spike deconvolution filter');
subplot(223),linesad(x_spike_deconv);title('x dec. with sp. dec. filter')
subplot(224),linesad(x_recovered);title('Convolution of de-spike x and wm')


[cc_obt_conv]=filter2coeff(peo_spike,lx-1);cc_obt_conv(1)=0;
cc_obt_conv=normalize(cc_obt_conv).*cc(it(1));

figure;
subplot(221);linesad(wav);
subplot(222);linesad(tt,R);title('Full Trace');
subplot(223);linesad(R_conv);title('Band limited trace');
subplot(224);linesad(tt,peo_spike);title('filter from band limited trace');
%unif_scale(gcf,222);
figure;
subplot(221);linesad(tt,cc);title('reflection coeff (given)');
subplot(222);linesad(tt,cc_obt);title('reflection coeff from Claerbout');
subplot(223);linesad(tt,cc_obt_conv);title('reflection coeff from band lim trace');
res=resta(cc,cc_obt_conv);
subplot(224);linesad(tt,res);title('Residuals r.c. :full trace - band lim. trace');
%unif_scale(gcf,221);

figure
for ii=1:3;
   delay=2*ii+14;
   [x_deconv_gap,peo_gap,w_gap,x_recovered_gap]=spikedecon(R_conv,bigdiag,delay);
   subplot(229+2*ii);linesad(peo_gap);title(strcat('peo-gap, delay=',int2str(delay)))
	subplot(230+2*ii),linesad(real(w_gap));title('w=peo-gap^-^1');
end;

delay=21;
[x_deconv_gap,peo_gap,w_gap,x_recovered_gap]=spikedecon(R_conv,bigdiag,delay);

figure,
subplot(221),linesad(peo_gap);title(strcat('peo-gap, delay=',int2str(delay)))
subplot(222),linesad(real(w_gap));title('w_gap=peo-gap^-^1');
subplot(223),linesad(x_deconv_gap);title('x deconvolved with gap deconvolution filter')
subplot(224),linesad(x_recovered_gap);title('x recovered with g0')

x_gap_dec=convlim(w_gap(1:delay),x_deconv_gap);

figure,
subplot(221),linesad(x_gap_dec);title(strcat('x-gap-dec, delay=',int2str(delay)))
subplot(222),linesad(wav(1:lw));title('wav')
subplot(224),linesad(x_deconv_gap(1:delay+lw));title('x deconvolved with gap deconvolution filter')
