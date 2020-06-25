function plothoriz(a,minf,df,mytitle,shift)
% plothoriz(a,minf,df,mytitle) 
% plot a matrix a horizontally wtih vertical axis increasing from
% top to bottom. 
% where
% a matrix to plot horizontally 
% minf min vertical axis
% df   vertical axis interval 
% title a string to put in the title 

[nr nc]=size(a);
a=a/max(max(abs(a))); %normalize

if (nargin<4) 
 mytitle='Figure';
 shift=1;
end

if (shift==1)
  for i=1:nc
%    [i*df i]
    a(:,i)=i*df+a(:,i);
  end
end


maxf=minf+nc*df;
figure,
if (shift==0)
  plot(a*maxf+minf),
else
  plot(a);
end

axis ij,
figure(gcf),

% Get handler to axes

% Change bottom axis to top
set(gca,'XAxisLocation','top');
% Delete left axis 
set(gca,'YTick',[0 10^6]);
% Change Font
set(gca,'FontSize',14);
% Limits
set(gca,'XLim',[0 nr]);
% Reverese X axis
set(gca,'XDir','reverse');

%axis([1 nr minf-df maxf+df])
xlabel('trace#')
ylabel('freq(Hz)')
title(mytitle)