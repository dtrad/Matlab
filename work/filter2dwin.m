function  [dz]=filter2dwin(d,f,F)
% Apply shaping filter in windows

% size of data
[nt nh]=size(d);

% size of filter
[nf1 nf2]=size(f);

% number of jumps in axis 1
nj1=ceil(nt/nf1);
ntz=nj1*nf1;

% number of jumps in axis 2
nj2=ceil(nh/nf2);
nhz=nj2*nf2;

% zero padding

dz=[d zeros(nt,nhz-nh)];
dz=[dz;zeros(ntz-nt,nhz)];

% Apply the filter in windows
for i=0:nj1-1
  for j=0:nj2-1
    w=window(dz,i*nf1+1,nf1,j*nf2+1,nf2);
    [w]=filter2d(w,f,F);
    dz=patching(w,dz,i*nf1+1,nf1,j*nf2+1,nf2);
  end
end


dz=dz(1:nt,1:nh);



