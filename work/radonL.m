function d=radonL(v,p,dt,h)

v=seis_shape(v);
[nt,np]=size(v);
nh=length(h);

if nargin<3|isempty(dt) dt=0.004;end 
if nargin<4|isempty(h) 
     nh=np;
     h=0:nh-1;
     h=h*25;
end
V=fft(v,nt); % data are zeropad to nt power of 2.
fs=1/dt;  % sample freq.
w=2*pi*(0:nt/2-1)*fs/nt; % angular frequency.
VH=V(1:nt/2,:);   
DR=zeros(1,nh);
nf=length(w);
for f=2:nf-30,
   [FH,F]=radonmat2(w(f),h,p);
   vh=VH(f,:).';  %'
   dr=F*vh;
   DR=[DR; dr.'];  %'
end

DR=zeropadm(DR,nt/2);
DD=duplic(DR);
d=ifft(DD);
t=0:nt-1;t=t*dt;
