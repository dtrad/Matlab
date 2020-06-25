% Program zero_offset
% Zero offset seismic reflection
% Compare programs Claerbout, Matlab (Levinson) and Porsani 
% GEOP 520: Direct Studies
% Daniel Trad: UBC. 
echo on
clear;
myclean;
dt=0.02;
lx=128;


t=[0.3 0.4];
rf=[0.5 -0.4];
n=length(rf);

it=round(t./dt+1);
cc=zeros(1,lx);cc(1)=1;
cc(it)=rf;
tt=0:dt:((lx-1)*dt);

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Claerbout's programs

[R,a]=rf2waves(cc,lx);
[cc_obt]=filter2coeff(a,lx-1);
a_lev=levinson(R,lx-1);cc(1)=0;  % PEO from Matlab

figure,
subplot(221);linesad(tt,R);title('waves, zero offset');
subplot(222);linesad(tt,a);title('filter');
subplot(223);linesad(tt,cc);title('reflection coeff (given)');
subplot(224);linesad(tt,cc_obt);title('reflection coeff (obtained)');

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Porsani's programs

[g,a_por]=levbak(lx,-cc);   
c_por=levfor(lx-1,a);c_por=-c_por;
R_por=acpeo(lx,lx,g);R_por=R_por(1:lx);

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Compare filter PEO from Claerbout and Porsani 

figure;
subplot(221);linesad(tt,a);title('filter from Claerbout');
subplot(222);linesad(tt,a_lev);title('filter from Levinson (peo))');
subplot(223);linesad(tt,a_por);title('filter from Porsani');
res=resta(a,a_por);
subplot(224);linesad(tt,res);title('Residuals from Claerbout-Porsani');
 
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Compare filter reflection coefficients from Claerbout and Porsani 
cc(1)=1;c_obt(1)=1;c_por(1)=1;
figure;
subplot(221);linesad(tt,cc);title('reflection coeff (given)');
subplot(222);linesad(tt,cc_obt);title('reflection coeff from Claerbout');
subplot(223);linesad(tt,c_por);title('reflection coeff from Porsani');
res=resta(cc,c_por);
subplot(224);linesad(tt,res);title('Residuals from Claerbout-Porsani');

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Compare trace from Claerbout and Porsani 

figure,
subplot(221),linesad(tt,R);title('R Claerbout')
subplot(222),linesad(tt,R_por(1:lx));title('R porsani')
subplot(223),linesad(resta(R,R_por));title('RC-R_por')

% Generate first and second order ringing multiples
cc_r=zeros(size(cc));cc_r(1)=1;cc_r(it(1))=rf(1);
[R_r,a_r]=rf2waves(cc_r,lx);

figure,
subplot(221);linesad(R_r);title('Rtemp_r')
subplot(222);linesad(a_r);title('a_r');
%subplot(223);linesad(R_sof);title('Rtemp_sof')
%subplot(224);linesad(a_sof);title('a_sof');

%temp1(nbig1+1)=xout(nbig1+1);
%temp2(nbig2+1)=xout(nbig2+1);
   
%p_m=convlim(a,temp1,lx,p_m);
%temp1=convlim(filt_1,temp2,lx,temp1);
%temp2=convlim(filt_1,temp1,lx,temp2);


% Compare multiple removal
R_noah=noahdecon(R);

figure,
subplot(221),linesad(tt,R);title('R Claerbout')
subplot(222),linesad(tt,R_noah(1:lx));title('R from noah deconvolution')
subplot(223),linesad(resta(R,R_noah));title('RC-R_noah')

fb=50;
lw=50;
wav=ricker(lw,fb,dt);
wn=-max(abs(wav));
wav=wav./wn;
R_temp=convlim(wav,R);

R_conv=R_temp(lw/2:lx+lw/2-1)';R_conv(1)=2.5;
a_lev_conv=levinson(R_conv,lx-1);
[cc_obt_conv]=filter2coeff(a_lev_conv,lx-1);

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


