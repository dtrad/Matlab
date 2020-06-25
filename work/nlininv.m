function [utr,vr,JD,JPF,J,sigmap]=nlininv(ur,UH,V,CI,epsilon,NITER,tolmodel,pnorm,C,axis1,optionc,initer,enditer,sigmap)
%	function [utr,vr,JD,JPF,J,sigmap]=lininv(ur,UH,V,CI,epsilon,NITER,tolmodel,pnorm,C,axis1,optionc,initer,enditer,sigmap)
% 	Daniel Trad May-27-98- UBC-EOS.
global w dte p h t np nh nt Power WV WU W dh h_near ff tor dp dt

for kk=initer:enditer;
      kk
      figure,
      
      if (initer~=1|enditer~=1) 
         sigmap(kk)=(10^(kk-1))*mean(sqrt(diag(C)));
      end   
      
      % Initial Model
      ii=1:np;
      tol(kk)=epsilon;
      if (optionc==1) V=zeros(size(V));end
      [V,JPF]=ForwNL(UH,V,tol(kk),NITER,tolmodel,pnorm,sigmap(kk),CI);
		V=V.';
		VD=duplic(V);
      vr=ifft(VD);
      [utr,JD(kk)]=backward(V,UH);
      if pnorm==10
         JP(kk)=sigmap(kk).*sum(sum(dp.*log(1+(abs(V(1:60,:))./sigmap(kk)).^2)));
      elseif pnorm==1
         JP(kk)=sigmap(kk).*sum(sum(dp.*(abs(V(1:60,:))./sigmap(kk)).^pnorm));
      else
         JP(kk)=sigmap(kk).^2.*sum(sum(dp.*(abs(V(1:60,:))./sigmap(kk)).^pnorm));
      end   
      J(kk)=JP(kk)+JD(kk);
      if Power==1 
         plotdata(ur,vr,utr,axis1)
      elseif Power==2
         plotdat2(ur,vr,utr,axis1)
      end   
		mytitle(sigmap(kk),pnorm);
end;
     	
plotnorm(sigmap,J,JP,JD);
mytitle2(pnorm);
%=====================================================================