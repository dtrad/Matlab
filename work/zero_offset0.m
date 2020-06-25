% Program zero_offset2
% Zero offset seismic reflection
% Stability of inversion.
% Generated with Claerbout's programs
% Theory Koehler and Taner Geophysics Vol42-No6 pag.1199-1206
% GEOP 520: Direct Studies
% Daniel Trad: UBC. 

echo on
clear;
myclean;
dt=0.004;
lx=256;
option=1;
bigdiag=2;
optionm='y'; % Option for multiple removal

optionw='no_w'
%optionw='minr'
%optionw='kais'
% Claerbout's program uses different sign than Robinson
% With rf=-rf we use Robinson's convention.
%t=[0.2 0.5];
%rf=[0.5 -0.3];rf=-rf; 
dt=1;
t=[25 30 50 75 87 100 105 108 125 137];
rf=[.4 -.1 .05 .10 -.03 -.15 .2 -.1 .15 .12];
n=length(rf);

it=round(t./dt+1);
cc=zeros(1,lx);
cc(it)=rf;
tt=0:dt:((lx-1)*dt);


[R,a]=rf2waves(cc,lx);
[cc_obt]=filter2coeff(a,lx-1);cc_obt(1)=0;
a_lev=levinson(R,lx-1);cc(1)=0;  % PEO from Matlab


% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

figure,
subplot(221);linesad(tt,R);title('waves, zero offset');
subplot(222);linesad(tt,a);title('filter');
subplot(223);linesad(tt,cc);title('reflection coeff (given)');
subplot(224);linesad(tt,cc_obt);title('reflection coeff (obtained)');
unif_scale(gcf,221);
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Compare filter PEO from Claerbout

figure;
subplot(221);linesad(tt,a);title('filter from Claerbout');
subplot(222);linesad(tt,a_lev);title('filter from Levinson (peo))');
% Remove multiples with noah deconvolution
unif_scale(gcf,221);
if optionm=='y' Rp=noahdecon(R);Rp(1)=0;end
%a_lev=levinson(Rp,lx-1);cc(1)=0;  % PEO from Matlab
%[cc_obt]=filter2coeff(a_lev,lx-1);cc_obt(1)=0;
%R=Rp;R(1)=1;
figure,
subplot(221);linesad(tt,Rp);title('primaries, zero offset ');
%subplot(222);linesad(tt,a);title('filter');
subplot(223);linesad(tt,cc);title('reflection coeff (given)');
%subplot(224);linesad(tt,cc_obt);title('reflection coeff (obtained)');
unif_scale(gcf,221);

alfa=1;
tau=0;
R=alfa*R;
R(1)=1;
if tau>=0 R=[R(1) zeros(1,tau) R(2:lx-tau)];
elseif tau<0 R=[R(1) R(-tau+2:lx)];R=padzeros(R,lx);
end;   

a_lev=levinson(R,lx-1);cc(1)=0;  % PEO from Matlab

[cc_obt]=filter2coeff(a_lev,lx-1);cc_obt(1)=0;

R=R./alfa;R(1)=1;
a_lev./alfa;a_lev(1)=1;
cc_obt=cc_obt./alfa;

figure,
subplot(221);linesad(R);title('waves, zero offset');
subplot(222);linesad(a_lev);title('filter');
subplot(223);linesad(cc);title('reflection coeff (given)');
subplot(224);linesad(cc_obt);title('reflection coeff (obtained)');
unif_scale(gcf,221);
if optionw~='no_w'
% Part II.
% Recovering the initial reflection coefficients
% from the trace after convolution with a wavelet.
% The following wavelets will be used:
% ricker
% ricker minimum phase
% kaiser

fb=50;
lw=30;
wav=seismwave(optionw,lw,dt,fb);lw=length(wav);
wav=normalize(wav);
R_temp=convlim(wav,R);

% Wrap around for NMP wavelet
if optionw=='rick' wa=11;else wa=1;end
R_conv=R_temp(wa:lx);
R_conv(lx-(wa-1):lx)=zeros(1,wa);

R_conv(1)=bigdiag;R_conv=R_conv./bigdiag; %stabilize the Levinson recursion 
a_lev_conv=levinson(R_conv,lx-1);a_lev_conv(lx)=0;R_conv(1)=0;
[cc_obt_conv]=filter2coeff(a_lev_conv,lx-1);cc_obt_conv(1)=0;
%cc_obt_conv=cc_obt_conv./max(abs(cc_obt_conv)).*cc(it(1));

R_conv=R_conv*bigdiag;R_conv(1)=1;
a_lev_conv=a_lev_conv*bigdiag;a_lev_conv(1)=1;
cc_obt_conv=cc_obt_conv*bigdiag;

figure;
subplot(221);linesad(tt,R_conv);title('Band limited trace- Non MP Ricker');
subplot(222);linesad(tt,a_lev_conv);title('filter from band limited trace');
subplot(223);linesad(tt,cc);title('reflection coeff (given)');
subplot(224);linesad(tt,cc_obt_conv);title('reflection coeff from band lim trace');
unif_scale(gcf,222);

end;

