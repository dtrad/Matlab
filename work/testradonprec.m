function [datartls,datarthr]=testradonprec(data,h,q,t,niter,method)

if (nargin<5) niter=3;end

nq=length(q);
dt=t(2)-t(1);
%method=1; % WTCGLS
vmax=[];
vtrue=[];
eps2=1e-1; 
step=0.95; 
curve='PRT';

[v,dr,q]=radon0prec(data,h,nq,q,vmax,dt,method,[],eps2,step,curve);
datartls=v;
%[v,dr,p]=radon0(d,h,np,vmin,vmax,dt,method,vtrue,eps2,step,rmethod)
figure,wigb(datartls)

for i=1:niter
  vtrue=v;
  [v,dr,q]=radon0prec(data,h,nq,q,vmax,dt,method,vtrue,eps2,step,curve);
end
datarthr=v;



