function unif_scale(n,g,ni,nf)
% unif_scale(n,g,ni,nf)
% Makes all four plots with uniform scale given by subplot(g) on figure(n)
% g must 221,222,223 or 224
% n maybe gcf or the number of the figure
% Daniel Trad 22/07/98
if nargin<4 ni=1;nf=4;end
figure(n)
subplot(g);v=axis;
for ii=ni:nf;subplot(220+ii);axis(v);end;