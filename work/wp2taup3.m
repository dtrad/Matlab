function [Rtc]=wp2taup2(R,wav,ppaxis,ttaxis,option,decim)
% wp2taup2(R,wav,ppaxis,ttaxis)
% R contains all frequencies if option='full'
% but the result is forced to be real.
% option='comp'; Uses all frequencies and output is complex;
% 
% Daniel Trad - UBC- 20-08-98
if nargin < 6 decim=0;end
if nargin>=5 option=option(1:4);end

if nargin < 5, option='half';end
R=seis_shape(R);
[NF,NP]=size(R);

if option=='half' NF=2*NF;end

if nargin < 4, ttaxis=1:NF;end
if nargin < 3, ppaxis=1:NP;end
if nargin < 2, wav=1;end

if option=='half' R=(R(1:NF/2,:));R=duplic(R);end
if option=='full' R=(R(1:NF/2,:));R=duplic(R);end

Rt=ifft(R);
if decim~=0
   Rtd=decimatem(Rt,decim);
   Rt=Rtd;clear Rtd;
	NF=NF/decim;
end
NF=length(ttaxis);

%for ii=1:NP,Rt(:,ii)=((0.95)^ii)*convlim(Rt(:,ii),wav,NF);end
for ii=1:NP,Rtc(1:NF,ii)=convlim(Rt(1:NF,ii),wav,NF);end
figure,wigb(Rtc,1,ppaxis,ttaxis);title('R from real freq');
xlabel('p (s/m)');ylabel('t (sec)');
