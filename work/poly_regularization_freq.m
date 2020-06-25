function poly_regularizatoin_freq(np,eps)
if (nargin < 2) eps=1e-3;end

load datairreg2.mat
hi=h;
[nt nh]=size(data);

dt=t(2)-t(1);

% generate h regular from h irregular
% (even if h is already regular
dh = (h(end)-h(1))/(length(h)-1);
hi=h;
h=hi(1):dh:hi(end);
if (length(h) ~= length(hi))
    display('error');
    return;
end
datar=zeros(size(data));
%zero padding and fft along time    
D=fft(data,2*nt);
fs=1/dt;nf=nt;
w=2*pi*(0:nf-1)*fs/nt;
DH=D(1:nf,:);     % data in f,x
MH=zeros(nf,nh);  % model in f,x
MHR=zeros(nf,nh); % model with regularization

for f=1:nf;
    slicer=real(DH(f,:));
    slicei=imag(DH(f,:));
    slicerp=mypoly_int(hi,slicer,h,np,eps);
    sliceip=mypoly_int(hi,slicei,h,np,eps);
    MH(f,:)=slicerp+i*sliceip;
end
M=duplic(MH);

model=ifft(M);model=model(1:nt,:);
figure;wigb(data,1,hi,t);
figure;wigb(model,1,h,t);



return