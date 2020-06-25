%%% testing of 1-D acoustic internal code


%%% build test data

%mod_data_3;


mod_data_4;


d=dw;

%%% calculate 1st order mult estimate

b3=Intcalc_b3(d,20,300);

 

tm=1:500;

clf;plot(d(tm));hold on;plot(b3(tm),'r');

d3=d+b3;

plot(d3(tm),'g');

%%%

 

%% Export to su files

write_subin('test_d',1,512,d(1:512))
write_subin('test_b3',1,512,b3(1:512))

%%% addheaders

suaddhead ns=512 <test_d |sushw key=dt a=2000 >test_d.su
suaddhead ns=512 <test_b3 |sushw key=dt a=2000 >test_b3.su

----------
X-Sun-Data-Type: default
X-Sun-Data-Description: default
X-Sun-Data-Name: Intcalc_b3.m
X-Sun-Charset: us-ascii
X-Sun-Content-Lines: 101

function b3=Intcalc_b3(d,del) 

% Calculates only a single term of Int WEMS

[temp,Nt]=size(d);

 
b3=zeros(size(d));
 

for m=1+del:Nt,
  
	test1 =zeros(size(d(1:m)));
	test1(1:m-del)=d(1:m-del);

	test2=fliplr(d(1:m));

	temp=conv(test1,test2);
 
	A2 =zeros(size(1:m));

	A2(1+del:m)=temp(m+del:2*m-1);
	
	temp=conv(d(1:m),A2);
 
	b3(m)=temp(m);

	

end



return  


% 	Old code 1997


	test1 =zeros(size(d));
	test1(1:m-del)=d(1:m-del);

	test2=fliplr(d);

	temp=conv(test1,test2);

	A2(1:Nt)= (temp(Nt:2*Nt-1));

	test1=zeros(size(1:Nt));
	test1(1+del:Nt)=A2(1+del:Nt);
 
	temp=conv(d,test1);

	b3(m)=temp(m);


%%%  Test version using maximum event separation Del

if Del>Nt 
	Del=Nt;
end

for m=1+del:Nt,
  
	test1 =zeros(size(d(1:m)));

	if (m<Del+1)  

		test1(1:m-del)=d(1:m-del);
	else
	
		test1(m-Del:m-del)=d(m-Del:m-del);
	end
	
	test2=fliplr(d(1:m));

%	temp is length 2*m-1 

	temp=conv(test1,test2);
 
%	cross correlation in temp centered at m;

	A2 =zeros(size(1:m));


	if (Del<m-1)

	  	A2(1+del:1+Del)=temp(m+del:m+Del);

 	else
		A2(1+del:m)=temp(m+del:2*m-1);

 	end	

	temp=conv(d(1:m),A2);
 
	b3(m)=temp(m);
	

end

----------
X-Sun-Data-Type: default
X-Sun-Data-Description: default
X-Sun-Data-Name: mod_data_3.m
X-Sun-Charset: us-ascii
X-Sun-Content-Lines: 95


%%%%%%%%%%%%%%%  Build data model to test acoustic internal demult

ipsilon=.1;
epsilon=.01;
apsilon=.001;
upsilon=.0001;
varupsilon=1e-15;

nt=512;
 
v1=2000;
v2=4000;
v3=2000;
v4=4000;
delta_t=25;

t_g=0;

dt=.002;
dz= (v1*dt);
 
 
% create w vector centered about DC

dw=2* pi /((nt-1)*dt);

wnyq = (nt+2)/2;

% sign is fixed to follow theoretical FT convention
w=-(-(wnyq-1):(wnyq-2))*dw;

R_1=(v2-v1)/(v2+v1);
R_2=(v3-v2)/(v3+v2);
R_3=(v4-v3)/(v4+v3);
T_1=(2*v2)/(v2+v1);
T_1p=(2*v1)/(v2+v1);
T_1=1-R_1^2;
T_2=T_1*(1-R_2^2);

Wfilt=fftshift(trap(1,2,80,100,nt)) ;

 
Z_1=exp(i*w*.5*delta_t*dt);

Z_2=exp(i*w*2.3*delta_t*dt);

%Z_1=exp(i*w*1*delta_t*dt);
%Z_2=exp(i*w*.8*delta_t*dt); 

Z_3=exp(i*w*2.7*delta_t*dt);

%D=Z.^2.*(R_1+T_1*T_1p*R_2*Z.^2-R_1*T_1*T_1p*R_2^2*Z.^4);
%R= Z_1.^2.*(R_1+T_1*R_2*Z_2.^2+T_2*R_3*Z_3.^2);

%R= Z_1.^2.*(R_1+T_1*R_2*Z_2.^2);

Z=Z_2;
R= Z_2.*(R_1+T_1*R_2*Z.^2 - T_1*R_1*R_2^2*Z.^4./(1+R_1*R_2*Z.^2));

d = real(ifft(fftshift( R)));

dw = real(ifft(fftshift( R.*Wfilt)));
 
 
 
t1=zeros(size(1:2*nt));
d_fs=t1;
t1(1:nt)=d;
d_fs(1:nt)=d;

for j=1:5;
	t1=-conv(t1,d);
	t1=t1(1:2*nt);
	d_fs=d_fs+t1;
end

d_fs=d_fs(1:nt);

d_fs(nt-10:nt)=0*d_fs(nt-10:nt);

temp=fftshift(fft(d_fs));

d_fsw=real(ifft(fftshift(temp.*Wfilt)));

 
d_fs=d_fs(1:nt);

%temp=conv(d_fs,filter);

%d_fsw=temp(1:nt);  




----------
X-Sun-Data-Type: default
X-Sun-Data-Description: default
X-Sun-Data-Name: mod_data_4.m
X-Sun-Charset: us-ascii
X-Sun-Content-Lines: 97


%%%%%%%%%%%%%%%  Build data model to test acoustic internal demult

ipsilon=.1;
epsilon=.01;
apsilon=.001;
upsilon=.0001;
varupsilon=1e-15;

nt=512;
 
v1=2000;
v2=4000;
v3=2000;
v4=4000;
delta_t=25;

t_g=0;

dt=.002;
dz= (v1*dt);
 
 
% create w vector centered about DC

dw=2* pi /((nt-1)*dt);

wnyq = (nt+2)/2;

% sign is fixed to follow theoretical FT convention
w=-(-(wnyq-1):(wnyq-2))*dw;

R_1=(v2-v1)/(v2+v1);
R_2=(v3-v2)/(v3+v2);
R_3=(v4-v3)/(v4+v3);
%T_1=(2*v2)/(v2+v1);
%_1p=(2*v1)/(v2+v1);
T_1=1-R_1^2;
T_2=1-R_2^2;

Wfilt=fftshift(trap(1,2,90,100,nt)) ;

delta_t=25;

Z_0=exp(i*w*1*delta_t*dt);
Z_1=exp(i*w*1.0*delta_t*dt);
Z_2=exp(i*w*3*delta_t*dt);

RR_2= (R_2+T_2*R_3*Z_2./(1+R_2*R_3*Z_2));

%RR_2= R_2 ;
%RR_2=0;

D= Z_0.*(R_1+T_1.*RR_2.*Z_1./(1+R_1.*RR_2.*Z_1));

%R= Z_0.*(R_1+T_1.*RR_2.*Z_1./(1+R_1.*RR_2.*Z_1));
 
%D= 1+ R./( 1+R ) ;

%D= 1+ R.*( 1-R +R.^2 -R.^3+R.^4) ;
 
temp = real(ifft(fftshift( D))); d=temp(1:nt);

d = real(ifft(fftshift( D)));
 
dw = real(ifft(fftshift( D.*Wfilt)));

 
t1=zeros(size(1:2*nt));
d_fs=t1;
t1(1:nt)=d;
d_fs(1:nt)=d;

for j=1:5;
	t1=-conv(t1,d);
	t1=t1(1:2*nt);
	d_fs=d_fs+t1;
end

d_fs=d_fs(1:nt);

d_fs(nt-10:nt)=0*d_fs(nt-10:nt);

temp=fftshift(fft(d_fs));

d_fsw=real(ifft(fftshift(temp.*Wfilt)));

 
d_fs=d_fs(1:nt);

%temp=conv(d_fs,filter);

%d_fsw=temp(1:nt);  




