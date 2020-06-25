% Program zero_offset2
% Zero offset seismic reflection- Wavelet Included.
% First and second order ringing.
% Generated with Claerbout's programs
% Theory Robinson and Treitel pag.288-289
% GEOP 520: Direct Studies
% Daniel Trad: UBC. 
echo on
clear;
myclean;
dt=0.004;
lx=256;
option=1;
%optionw=='rick'
optionw='minr'
%optionw='kais'

% Claerbout's program uses different sign than Robinson
% With rf=-rf we use Robinson's convention.
t=[0.2 0.5];
rf=[0.5 -0.3];rf=rf; 
n=length(rf);

it=round(t./dt+1);
cc=zeros(1,lx);
cc(it)=rf;
tt=0:dt:((lx-1)*dt);

% First we compute the first and second order ringing filter
cctemp=zeros(size(cc));Rtemp=zeros(lx,1);
cctemp(it(1))=cc(it(1));
[Impulse,filter_1st]=rf2waves(cctemp,lx);
filter_2nd=convlim(filter_1st,filter_1st,lx);
Impulse2=convlim(Impulse,Impulse,lx);

[R,a]=rf2waves(cc,lx);
[cc_obt]=filter2coeff(a,lx-1);cc_obt(1)=0;
a_lev=levinson(R,lx-1);cc(1)=0;  % PEO from Matlab

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

figure
%subplot(221);linesad(tt,Impulse);title('Impulse Response for first order reverberation at zero offset');
%subplot(222);linesad(tt,Impulse2);title('Impulse Response for second order reverberation at zero offset');
%subplot(223);linesad(tt,filter_1st);title('first order ringing filter (1/Impulse Response)');
%subplot(224);linesad(tt,filter_2nd);title('second order ringing filter (1/Impulse Response)^2');

subplot(221);linesad(tt,Impulse);title('(a)');
subplot(222);linesad(tt,Impulse2);title('(b)');
subplot(223);linesad(tt,filter_1st);title('(c)');
subplot(224);linesad(tt,filter_2nd);title('(d)');


figure,
%subplot(221);linesad(tt,R);title('waves, zero offset');
%subplot(222);linesad(tt,a);title('filter');
%subplot(223);linesad(tt,cc);title('reflection coeff (given)');
%subplot(224);linesad(tt,cc_obt);title('reflection coeff (obtained)');

subplot(221);linesad(tt,R);title('(a)');
subplot(222);linesad(tt,a);title('(b)');
subplot(223);linesad(tt,cc);title('(c)');
subplot(224);linesad(tt,cc_obt);title('(d)');


% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Compare filter PEO from Claerbout and Porsani 

figure;
subplot(221);linesad(tt,a);title('filter from Claerbout');
subplot(222);linesad(tt,a_lev);title('filter from Levinson (peo))');
 
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Filtering  
% Matlab function filter works as 
%    Y = FILTER(B,A,X) filters the data in vector X with the
%    filter described by vectors A and B to create the filtered
%    data Y.  The filter is a "Direct Form II Transposed"
%    implementation of the standard difference equation:
% a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                          - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% For first order filter
% The impulse response is 1/(1+c1.z^tau1). We want to deconvolve this
% which is equivalent to filter by (1+c1.z^tau1).
% We have R=x/(1+c1.z^tau1) so that x=R (1+c1.z^tau1)

R(1)=0;
R_1st=filter(filter_1st,1,R);

% For second order filter:
if option==1 
R_1st(it(1))=0; % eliminate the primary
R_2nd=filter(filter_1st,1,R_1st);
R_2nd(it(1))=R(it(1)); % add back the primary
elseif option==2
% Othe way is
R_1st(it(1))=0; % eliminate the primary
R_2nd=filter(filter_2nd,1,R_1st);
R_2nd(it(1))=R(it(1)); % add back the primary
end;
R_1st(it(1))=R(it(1)); % add back the primary

figure,
%subplot(221);linesad(R);title('R')
%subplot(222);linesad(R_1st);title('R_1st');
%subplot(223);linesad(R_2nd);title('R_2nd')
%subplot(224);linesad(resta(R,R_2nd));title('Residuals');

subplot(221);linesad(R);title('(a)')
subplot(222);linesad(R_1st);title('(b)');
subplot(223);linesad(R_2nd);title('(c)')
subplot(224);linesad(resta(R,R_2nd));title('(d)');

% Compare multiple removal
R_noah=noahdecon(R);

figure,
%subplot(221),linesad(tt,R);title('R Claerbout')
%subplot(222),linesad(tt,R_noah(1:lx));title('R from noah deconvolution')
%subplot(223),linesad(resta(R,R_noah));title('RC-R_noah')

subplot(221),linesad(tt,R);title('(a)')
subplot(222),linesad(tt,R_noah(1:lx));title('(b)')
subplot(223),linesad(resta(R,R_noah));title('(c)')

% Part II.
% Recovering the initial reflection coefficients
% from the trace after convolution with a wavelet.
% The following wavelets will be used:
% ricker
% ricker minimum phase
% kaiser

fb=25;
lw=50;
if optionw=='rick'
	wav=ricker(lw,fb,dt);
	figure,plot(wav)
elseif optionw=='minr'
   load minricke.dat;lw=length(minricke);
   wav=minricke(6:lw-5);
   figure,plot(wav)
elseif optionw=='kais'
   load kaiserspike;lw=length(temp);
   inipoint=50;wav=temp(inipoint:lw-inipoint+1);
   %wav=fftshift(temp(1:lw/2));
   figure,plot(wav)
end

wn=max(abs(wav));
wav=wav./wn;
R(1)=1;R_temp=convlim(wav,R);
primary1=zeros(lx,1);primary1(it(1))=cc(it(1));
primary1=convlim(wav,primary1,lx);
ref=zeros(lx,1);ref(1)=1;
ref=convlim(wav,ref,lx);
%R_conv=R_temp(lw/2:lx+lw/2-1)';
R_conv=R_temp(1:lx);

figure,
subplot(221);plot(primary1);title('primary1');
subplot(222);plot(ref);title('ref');
subplot(223);plot(R_conv);title('R_conv');

%R_conv(1)=1; %stabilize the Levinson recursion 
a_lev_conv=levinson(R_conv,lx-1);R_conv=R_conv-ref;
[cc_obt_conv]=filter2coeff(a_lev_conv,lx-1);cc_obt_conv(1)=0;
cc_obt_conv=cc_obt_conv./max(abs(cc_obt_conv)).*cc(it(1));

figure;
subplot(221);linesad(wav);
subplot(222);linesad(tt,R);title('Full Trace');
subplot(223);linesad(R_conv);title('Band limited trace');
subplot(224);linesad(tt,a_lev_conv);title('filter from band limited trace');

figure;
subplot(221);linesad(tt,cc);title('reflection coeff (given)');
subplot(222);linesad(tt,cc_obt);title('reflection coeff from Claerbout');
subplot(223);linesad(tt,cc_obt_conv);title('reflection coeff from band lim trace');
res=resta(cc,cc_obt_conv);
subplot(224);linesad(tt,res);title('Residuals r.c. :full trace - band lim. trace');

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Filtering the trace convolved with wavelet
% For first order filter


R_conv_1st=filter(filter_1st,1,R_conv);

% For second order filter:
if option==1 
R_conv_1st=R_conv_1st-primary1; % eliminate the primary
R_conv_2nd=filter(filter_1st,1,R_conv_1st);
R_conv_2nd=R_conv_2nd+primary1; % add back the primary
elseif option==2
% Othe way is
R_conv_1st=R_conv_1st-primary1; % eliminate the primary
R_conv_2nd=filter(filter_2nd,1,R_conv_1st);
R_conv_2nd=R_conv_2nd+primary1; % add back the primary
end;
R_conv_1st=R_conv_1st+primary1; % add back the primary

figure,
subplot(221);linesad(R_conv);title('R')
subplot(222);linesad(R_conv_1st);title('R_1st');
subplot(223);linesad(R_conv_2nd);title('R_2nd')
subplot(224);linesad(resta(R_conv,R_conv_2nd));title('Residuals');

R_decon=deconv(R_conv,wav);
R_decon_1st=deconv(R_conv_1st,wav);
R_decon_2nd=deconv(R_conv_2nd,wav);

figure,
subplot(221);linesad(R_decon);title('R_decon')
subplot(222);linesad(R_decon_1st);title('R_decon_1st');
subplot(223);linesad(R_decon_2nd);title('R_decon_2nd')
subplot(224);linesad(resta(R_decon,R_decon_2nd));title('Residuals');

% Compare multiple removal
R_conv_noah=noahdecon(R_conv,wav);

figure,
subplot(221),linesad(tt,R_conv);title('R Claerbout')
subplot(222),linesad(tt,R_conv_noah(1:lx));title('R from noah deconvolution')
subplot(223),linesad(resta(R_conv,R_conv_noah));title('RC-R_noah')

R_decon_noah=deconv(R_conv_noah,wav);

figure,
subplot(221);linesad(R_decon);title('R_decon')
subplot(222);linesad(R_decon_noah);title('R_decon_noah');
subplot(224);linesad(resta(R_decon,R_decon_noah));title('Residuals');

 