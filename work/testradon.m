function [datartls,datarthr]=testradon(data,h,q,t,niter)

if (nargin<5) niter=3;end

nq=length(q);
dt=t(2)-t(1);
method=3; % WTCGLS
vmax=[];
vtrue=[];
eps2=1e-3; 
step=0.95; 
curve='PRT';

[v,dr,q]=radon0b(data,h,nq,q,vmax,dt,method,[],eps2,step,curve);
datartls=v;
%[v,dr,p]=radon0(d,h,np,vmin,vmax,dt,method,vtrue,eps2,step,rmethod)
for i=1:niter
  vtrue=v;
  [v,dr,q]=radon0b(data,h,nq,q,vmax,dt,method,vtrue,eps2,step,curve);
end
datarthr=v;



