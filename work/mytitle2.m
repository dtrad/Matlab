function mytitle2(pnorm)
if pnorm==1 TXT=sprintf('Norm=L1');end
if pnorm==2 TXT=sprintf('Norm=L2');end
if pnorm==10 TXT=sprintf('Norm=Cauchy');end
subplot(221);text(1,700,TXT)