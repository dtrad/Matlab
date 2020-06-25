function poly_regularization(np,eps)
if (nargin < 2) eps=1e-3;end

load datairreg.mat
data=data(:,1:20);
h=h(1:20);
hi=h;
[nt nh]=size(data);
amax=max(max(data));
data=data+amax*0.1*rand(size(data));

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
figure;wigb(data,1,hi,t);

for it=1:nt
    slice=data(it,:);
    slicep=mypoly_int(hi,slice,h,np,eps);
    datar(it,:)=slicep;
end


figure;wigb(datar,1,h,t);



return