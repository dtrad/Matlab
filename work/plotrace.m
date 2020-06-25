function plotrace(data,dh,nh,h_near,t,dt)
%function plotrace(data,dh,nh,h_near,t)
% Daniel Trad - May 8/98  UBC-EOS
if nargin<6 dt=0.004;end;
if nargin<5 [nt nc]=size(data); t=(0:(nt-1))*dt;end;
if nargin<4, h_near=0;end
if nargin<3, nh=nc;end
if nargin<2, dh=50;end

data=real(data).*dh;
for jj=1:nh, plot(t,data(:,jj)+(jj-1)*dh+h_near),hold on,end,hold off
