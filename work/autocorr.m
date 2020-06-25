function [Rx]=autocorr(in,option,dh,h_near)
% function [Rx]=autocorr(in,option,dh,h_near)
% input:
% in (trace)
% option; plot 'y' or 'n' (default 'n')
% dh: offset interval (default 50)
% h_near: (defaul=0;)
% output: autocorrelation for every trace;
%
% Daniel Trad-- June 2, 1998- UBC.

if nargin<2, option='n';end
if nargin<3, dh=50;end
if nargin<4, h_near=0;end

[nrow ncol]=size(in);
if nrow<ncol in=in.';end
[nrow ncol]=size(in);

for hh=1:ncol;
   Rx(:,hh)=xcorr(in(:,hh),'coeff');
end;   

t=1:nrow;

if option=='y'
   plotrace(Rx(nrow:2*nrow-1,:),dh,ncol,h_near,t)
   xlabel('lag');ylabel('offset');title('Autocorrelation function')
   figure(gcf)
end
