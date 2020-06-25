function mytitle(sigmap,pnorm)
if pnorm==1 TXT=sprintf('Sigmap=%15.2e Norm=L1',sigmap);end
if pnorm==2 TXT=sprintf('Sigmap=%15.2e Norm=L2',sigmap);end
if pnorm==10 TXT=sprintf('Sigmac=%15.2e Norm=Cauchy',sigmap);end
subplot(221);text(1,700,TXT)