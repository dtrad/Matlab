function  [dz]=ridgefilterwin(d,n1,n2,perc)
% Apply ridge filter in windows of size n1 x n2

% size of data
[nt nh]=size(d);

% number of jumps in axis 1
nj1=ceil(nt/n1);
ntz=nj1*n1;

% number of jumps in axis 2
nj2=ceil(nh/n2);
nhz=nj2*n2;

% zero padding

dz=[d zeros(nt,nhz-nh)];
dz=[dz;zeros(ntz-nt,nhz)];

% Apply the filter in windows
for i=0:nj1-1
  for j=0:nj2-1
    w=window(dz,i*n1+1,n1,j*n2+1,n2);
    [w]=ridgefilter(w,perc);
    dz=patching(w,dz,i*n1+1,n1,j*n2+1,n2);
  end
end


dz=dz(1:nt,1:nh);



