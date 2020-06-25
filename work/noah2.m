
%	program Noah

% ****** noah is an update to pandm to incorporate Noah
% ****** the removal of all surface multiples is the quest.
% ****** I use the f domain
% ------------ steps to compute reflection coeffs from seismogram
% ****** x is the seismogram
% ****** ppg are the refl. coeffs
%	   CALL PEO(NCF,G(2),X,POT) -------- computes the peo g
%	   g(1)=1
%	   call LEVfor(lcf,G,gpp)   -------- computed the refl. coeffs.
%   layerf obtains the synthetic seismogram 1+2R(z) from
%   the reflection coefficients which are read from 2
%   the output is read into 1
clear;
for i=1:20;close;end;
tadcomp='n';

if tadcomp=='y'
load g1tad;
load xtad;
load gpptad;
load xnoah;

xtad=x;
g1tad=g1;
gpptad=gpp;
xnoahtad=xnoah;

clear x;
clear g1;
clear gpp;
clear xnoah;
end;

dt=0.004;
optionfof='y';
optionsof='y';
iwav=0;
lx=300;
lxmax=500;
x=zeros(lxmax,1);
g=zeros(lxmax,1);
gpp=zeros(lxmax,1);
g2=zeros(lxmax,1);
x=zeros(lxmax,1);
xs=zeros(lxmax,1);
xin=zeros(lxmax,1);
xout=zeros(lxmax,1);


tt=[.2 .7];  
gt=[.4 -.3];

n=length(tt);
pp(1:n)=gt(1:n);

cc=gt(1);
nn=fix(tt(1)/2/dt+0.5);

% now set up the layered system


% first determine one-way time, then determine the step in delt of 4ms.

itt(1)=0
for i=1:n
	itt(i+1)=fix(tt(i)/2/dt+0.5);
	ktt=itt(i+1);
	if(i==1)nbig1=ktt;end
	gpp(ktt)=gt(i);
	g(ktt)=-gpp(ktt);
	pp(ktt)=g(ktt);
	for k=itt(i)+1:itt(i+1)-1
		gpp(k)=0.0;
		g(k)=gpp(i);
		pp(k)=g(k);
	end
end

figure,
subplot(221);plot(gpp),title('gpp');
subplot(222);plot(g),title('g');
subplot(223);plot(pp);title('pp');

nkeep=n;
n=ktt;
nbig2=n;

display('BIG n =');n
display('nbig1 =');nbig1
display('nbig2 =');nbig2

nfft=128;

cxin=zeros(nfft,1);
cxnoah=zeros(nfft,1);
ctemp=zeros(nfft,1);

% ****** compute the times for a check
ipos=0;
kk=0;

for i=1:n
	  ipos=ipos+1;
	  if(abs(g(i))~=0)
	  		kk=kk+1;
	  		tt(kk)=ipos*2.*dt;
	  		pp(i)=g(i);
	  end
end

xp(1)=1;

for i=2:nkeep+1
	ipos=itt(i)+1;
	xp(ipos)=-gt(i-1);
end

xp(1)=-1;
xp=-xp;

for ipass=1:2;

% ****** first determine the reverberation 1st order filter
	if(ipass==1)n=nbig1;end
	if(ipass==2)n=nbig2;end
   
	
	ncf=lx;
	%c=resta(gpp,gpptad,'y','gpp-gpptad');
	g(1)=1;
   [g1,g(2:n+1)]=levbak(n,gpp);
   %Check levfor
   rcoef=levfor(n,g);
   subplot(224);plot(rcoef);title(strcat('rcoef from levfor, ipass=',int2str(ipass)))
   figure
   if ipass==2
      subplot(221),plot(g1,'b');%hold on,plot(g1tad,'r');
      title('g1');
   subplot(222),plot(g);title('g');
   end   
   x=acpeo(n,lx,g1);
   
   if ipass==2
      aa=levinson(x,200);
      figure,
      subplot(221),plot(x,'b');title('output from acpeo');
      gpptemp=[1;gpp(:)];rxx_gpp=autocorr(gpptemp);lgpp=length(gpptemp);
      subplot(222),plot(rxx_gpp(lgpp:-1:lgpp./2+1));title('autocorr of gpp');
      subplot(223),plot(aa);title('peo')
      
   end;   
   %c=resta(g1,g1tad,'y','g1-g1tad');
   %c=resta(x,xtad,'y','x-xtad');
   
	lx=length(x);i=1:lx;
	xs(i)=x(i);
	pp(i)=x(i);

	xs(1)=-1;
	xs(i)=-xs(i);
	xin(i)=xs(i);
	xout(i)=xs(i);
	
	if(ipass==1)
		i=1:lx;
		filt_1(i)=xout(i+nn);
		xn=max(abs(filt_1));
   	filt_1=filt_1./xn;
	end
end

% ****** now produce the two primary reverberation series
	p_m(i)=zeros(lx,1);
	temp1(i)=zeros(lx,1);
	temp2(i)=zeros(lx,1);

% ****** for 1st order ringing and 2nd order ringing

	temp1(nbig1+1)=xout(nbig1+1);
   temp2(nbig2+1)=xout(nbig2+1);
   
   p_m=convlim(filt_1,temp1,lx,p_m);
   temp1=convlim(filt_1,temp2,lx,temp1);
   temp2=convlim(filt_1,temp1,lx,temp2);
   	
	temp2(1)=1;
   figure,
   subplot(221);plot(p_m);title('p_m');
   subplot(222);plot(temp1);title('temp_1');
   subplot(223);plot(temp2);title('temp_2');vv=axis;
   subplot(224);plot(filt_1);title('filt_1');
   subplot(221);axis(vv);
   subplot(222);axis(vv);
   subplot(224);axis(vv);
   
% ****** now sum

	p_m=p_m+temp2;
	xs=xin;
	
	xout(1)=1;

% ****** p_m is the 1st order and 2nd order series + primaries
% ****** now compute the Noah deconvolved series
% ****** the series is 1+R(z) and is obtained as
%	 1+R(z)=1+(X(z)-1)/(2-X(z))

	cxin=fft(xin);nfft=length(cxin);

	i=1:nfft;
   cxnoah(i)=1+(cxin(i)-1)./(2.-cxin(i));
   cxnoah=duplic(cxnoah(1:nfft./2));
   
   xnoah=real(ifft(cxnoah));
   %c=resta(xnoah,xnoahtad,'y','xnoah-xnoahtad');
   
   figure,
   subplot(221);plot(p_m);title('1st and 2nd order series');
	subplot(222);plot(xnoah);title('Noah deconvolved series');
	aa=levinson(xnoah,200);
	subplot(223),plot(aa);title('peo')

% ****** now compare the traces
   dref=zeros(lx,1); 
	dref(1)=1;
   xs=xs(1:lx);
   if iwav==1 
   	fb=25;
		dt=0.004;
   	lw=50;
   	wav=ricker(lw,fb,dt);
   	wn=max(abs(wav));
   	wav=wav./wn;
      xss=convlim(wav,xs,lx);
   	xs=xss;
      xin=xss;
      pp0=convlim(wav,xp,lx);
      dwav=convlim(wav,dref,lx);
      dref=dwav;
      subplot(223);plot(wav);title('wavelet');
      i=1:lx;
  		xp(i)=pp0(i);
   end
   
   
   subplot(224);plot(xs);title('xs');
      

	gpp=-gpp;
	gt=-gt;
   
   if(optionfof=='y')
	 	nn1=nn;
	 	cc1=cc;
		xss=xs-dref;
		% ------ now do a 1st order recursive filter.
		gpp(nn1)=cc1;

		for i=1:lx
			if(i-nn1<1)
				yy1=0;
			else
				yy1=xss(i-nn1);
			end
			pp0(i)=xss(i)+gpp(nn1)*yy1;
		end
      pp0=pp0(:);
		xout=pp0;
		xout=xout+dref;
            
		figure,
		subplot(221);plot(xs);title('Seismic trace');vv=axis;
		subplot(222);plot(xout);title('After 1st order');axis(vv);

		figure(gcf+1);subplot(221);plot(gpp);title('Gpp 1nd order');

	end

if(optionsof=='y')
   
%   define the first PRIMARY to be subtracted from i/p
%	 and call it primary1


		primary1=zeros(lx,1);
			
		if(iwav==1)lw=50;end
		if(iwav==0)lw=1;end
		for i=lw/2+nn1-lw/2:lw/2+nn1+lw/2
         if i>0 primary1(i)=xout(i)-dref(i);end
		end
		xss=xout-dref-primary1;
      subplot(223);plot(xout);title('xss');axis(vv);
      subplot(224);plot(primary1);title('primary');axis(vv);
      
% ------ now do a 2nd order recursive filter.
		gpp(nn1)=cc1;

		for i=1:lx
			if(i-nn1<1)
				yy1=0;
			else
				yy1=xss(i-nn1);
      	end
         pp0(i)=xss(i)+gpp(nn1)*yy1;
      end

		xout=pp0;
      xout=xout+dref+primary1;
      figure(gcf-1);      
      subplot(223);plot(xout);title('Seismic trace after 2nd order');axis(vv);
      subplot(224);plot(resta(xs,xout));title('residuals');axis(vv);
      figure(gcf+1);
      subplot(222);plot(gpp);title('Gpp 2nd order');
		aa=levinson(xout,200);
		subplot(223),plot(aa);title('peo')

end










