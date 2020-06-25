function [y,peo,w,xx]=spikedecon(x,bigdiag,delay,nn,option)
% Computes the PEO for a trace X
% the predicted spike or gap trace Y
% and the wavelet W=1/PEO
%		[y,peo]=spikedecon(x,bigdiag)
% Theory: Tad Ulrych: Notes from GEOP520b- UBC-CA
% x(t) trace; w(t) wavelet; r(t) reflectivity series (assumed white)
% x(t)=w(t)*r(t) (* stands for convolution)
% X(z)=W(z).R(Z)
% R(z).R(1/z)=sigma^2;
%
% X(z).X(1/z)=Rxx(z)=sigma^2.W(z).W(1/z)
% W(z)=1/G(z) with G(z)=PEO
% Rxx(z)=sigma^2.1/G(z).1/G(1/z)
%------------------------------
% Rxx(z).G(z)=sigma^2.1/G(1/z) |
%------------------------------
% and X(z).G(z)=R(z)--> reflectivity series!!!!
% ----------------------------------------------
% Daniel Trad- UBC- 30-07-98

lx=length(x);

if nargin<5 option='n';end
if nargin<4 nn=lx;end
if nargin<3 delay=1;end
if nargin<2 bigdiag=1;end  % Increases diagonal in Toeplitz matrix Rxx

ref=zeros(1,2*nn);ref(1)=1;

Rxx=autocorr(x(1:lx)); 
Rxx=Rxx(lx:2*lx-1);
Rxx(1)=bigdiag;

% Rxx(z).G(z)=[sigma^2,0,0,..] solved by Levinson recursion
peo=levinso1(Rxx,nn-delay,delay);

% If delay~=1 (gap pred.) zeros are set in the peo
la=length(peo);peo(la:la)=0;
az=zeros(1,delay-1);
peo=[peo(1) az -peo(2:la)];

w=deconv(ref,peo);  % Polinomial division
%w=deconv_freq_domain(ref,peo);  % Polinomial division
y=convlim(x,peo,lx); % R(z)=X(z).G(z)
xx=convlim(y,w,lx);  % X(z)=W(z).R(z)

if option=='y'
figure,
%subplot(221),linesad(Rxx);title('(a) Autocorrelation');
%subplot(222),linesad(peo);title('(b) spike deconvolution filter');
%subplot(223),linesad(x);title('(c) input trace')
%subplot(224),linesad(y);title('(d) output trace')

subplot(221),linesad(Rxx);title('(a)');
subplot(222),linesad(peo);title('(b)');
subplot(223),linesad(x);title('(c)')
subplot(224),linesad(y);title('(d)')

end