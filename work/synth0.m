function synth0(Stest)
% model
%close all
options='y'; %save
options2='n'; %save to binary file (float data)
optionm='n'; %multiples
%optionw='minr'; % wavelet
optionw='rick';
%optionw='freq';
optionn='n'; % noise;
offset='regu';
%offset='ireg';
%curve='line';
%curve='para';
curve='hype';
curve='mult';
curve='shif'

scalep=10;

pathfile='';
file='model5';
filename=[pathfile file];

tf=6.0;
dt=0.004;
dh=50;
nh=64;
nhfin=64;
h_near=0;
nw=200;
f=100;
aten=-1.2;
aten=0;
ngap=0;

S=(ones(5,1));
v(1)=1500;t0(1)=1.0;cr(1)=1;S(1)=1;
v(2)=2500;t0(2)=0.8;cr(2)=(1-cr(1)^2)*0.5;cr(2)=0;
v(3)=1500;t0(3)=1.0;cr(3)=-1;S(3)=Stest;
v(4)=4000;t0(4)=1.2;cr(4)=0.;
v(5)=4200;t0(5)=1.4;cr(5)=0.;


if (length(v)~=length(cr)) 
   display('Check the v and cr lenght');
end

% multiples
if optionm=='y'
v(5)=v(1);t0(5)=2*t0(1);cr(5)=-cr(1)^2;
v(6)=v(1);t0(6)=3*t0(1);cr(6)=+cr(1)^3;
v(7)=v(1);t0(7)=4*t0(1);cr(7)=-cr(1)^4;
v(8)=v(1);t0(8)=5*t0(1);cr(8)=+cr(1)^5;
v(9)=v(2);t0(9)=2*t0(2)+t0(1);cr(9)=-cr(2).^2*cr(1);
v(10)=v(2);t0(10)=t0(2)+2*t0(1);cr(10)=-cr(2)*cr(1).^2;
%v(11)=v(2);t0(11)=t0(3)+3*t0(1);cr(11)=-cr(3)*cr(1).^3;
%v(12)=v(1);t0(12)=1.8;cr(12)=+cr(1)^9;
%v(13)=v(1);t0(13)=2.0;cr(13)=+cr(1)^10;
end
nr=max(size(v));

if optionw=='rick'
	w=ricker(f,dt);%w=padzeros(w,nw);
	%figure,plot(w)
elseif optionw=='minr'
  	 load minricke.dat;
	w=minricke;
elseif optionw=='freq'
	tw=0:nw-1;tw=tw*dt;
	w=sin(2*pi*tw*f);

	wind=hanning(nw/4);
	wind=wind(:).';        %'
	wind=[wind(1,1:nw/8) ones(1,3*nw/4) wind(1,nw/8+1:nw/4)];
	w=wind(:).*w(:);
	%figure,plot(w);   
end


nt=round(tf/dt)+1;
xx=zeros(2*nt,2*nh);
if (offset=='ireg')
	h=rand_offset(nh,dh,scalep);h=h-h_near;
elseif (offset=='regu')
	h=0:nh-1;h=h*dh;h=h-h_near;
end




for rr=1:nr;
   for hh=1:nh
      %h(hh)=dh*(hh-1)+h_near;
      avo=abs(cos((h-h_near)./((nh-1)*dh)*pi));
      
      if (curve=='hype')
        	t=sqrt(t0(rr).^2+(h(hh)/v(rr)).^2);
      elseif (curve=='para')	
     	 	t=t0(rr)+(h(hh)/v(rr)).^2;
      elseif (curve=='line')
			t=t0(rr)+h(hh)/v(rr);
      elseif (curve=='mult')                 
            t=multif(t0(rr),h(hh),v(rr),0);
      elseif (curve=='shif')
          t=((S(rr)-1)/S(rr))*t0(rr)+sqrt((t0(rr)/S(rr)).^2+(((h(hh)/v(rr)).^2)/S(rr)));
      end;
	
      tt=(t./dt)+1;
      tt1=tt-floor(tt);
      tt2=1-tt1;
   	%if hh==23 tt, end
      %xx(tt,hh)=(xx(tt,hh)+cr(rr)).*avo;
     	xx(floor(tt),hh)=(xx(floor(tt),hh)+tt2*cr(rr));
     	xx(ceil(tt),hh)=(xx(ceil(tt),hh)+tt1*cr(rr));
	end;
end;
tt=1:nt;
w=w(:).';     %'
for hh=1:nh
      xxtemp=conv(xx(tt,hh),w);
      xxc(1:nt,hh)=xxtemp(1:nt);
      %xxc(1:nt,hh)=xxtemp(fix(nw/2)+1:nt+fix(nw/2));
end;

if optionn=='y';
   noise=randn(size(xxc));
   xxc=xxc+(noise-0.5)*0.003;
end;

xxtemp=xx./max(max(abs(xx)));
t=(0:nt-1)*dt;ntlim=ceil(nt/2);
figure,wigb(real(xx(1:ntlim,1:nh)),1,h,t(1:ntlim));
%figure;plotrace(xxtemp,dh,nh,h_near)
figure,wigb(real(xxc(1:ntlim,1:nh)),1,h,t(1:ntlim));

xxtemp=xxc./max(max(abs(xxc)));
%figure;plotrace(xxtemp(1:nt,:),dh,nh,h_near)
ylabel('time (sec)');xlabel('offset (m)');
title('Synthetic Seismic CDP');
if options=='y' x=xxc(1:256,1:40);h2=h(1:40);t2=t(1:256); save xxc3.mat x h2 dt t2;end;
if options2=='y';
   naux=(nh-1)/2;
   %xxz=[xxc(1:512,1:(naux-ngap+1)) xxc(1:512,(naux+1)-ngap/2) xxc(1:512,(naux+1)) xxc(1:512,(naux+1+ngap/2)) xxc(1:512,(naux+ngap:nh))];
   %xxz=[xxc(1:512,1:(naux+15)) xxc(1:512,(naux+25:nh))];
   %h=[h(1:(naux-ngap+1)) h(naux+1-ngap/2) h(naux+1) h(naux+1+ngap/2) h((naux+%xxz=[xxc(1:512,1:(naux+15)) xxc(1:512,(naux+25:nh))];
   %h=[h(1:(naux-ngap+1)) h(naux+1-ngap/2) h(naux+1) h(naux+1+ngap/2) h((naux+ngap:nh))];
   xxz=xxc(1:512,1:nhfin);
   h=h(1:nhfin);

   figure,wigb(xxz,1,h,t(1:512));   
   xxz=xxz(:);
   h=h(:);
   
   fid = fopen([filename '.bin'],'wb')
   for it=1:length(xxz)
         fwrite(fid,xxz(it),'float32');
   end
   fclose(fid)
   
   fid = fopen([filename '.off'],'wb')
   for it=1:length(h)
      fwrite(fid,h(it),'float32');
   end
   fclose(fid)
end;

return






