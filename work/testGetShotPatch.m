nfig=1;
for N=1:100:1000
    [rpatchX,rpatchY,thisShotX,thisShotY,shotFold]=getShotPatch(shotUniqueId,shotBinX,shotBinY,rcvrBinX,rcvrBinY,shotUniqueId(N));
    shotFold; 
    figure(nfig)
    plot(rpatchX,rpatchY,'o',thisShotX,thisShotY,'+');text= sprintf('shot patch ID = %d with fold = %d',shotUniqueId(N),shotFold);title(text); 
    nfig= nfig+1
end