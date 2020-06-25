function [utr,vr,V,JD,JP,J,sigmap]=lininv(ur,UH,CI,sigmap,pnorm,C,axis1,initer,enditer)

%  function [utr,vr,V,JD,JP,J,sigmap]=lininv(ur,UH,CI,sigmap,pnorm,C,axis1,initer,enditer)
% 	Daniel Trad May-27-98- UBC-EOS.

global w dte p h t np nh nt Power WV WU W dh h_near ff tor dp dt
for kk=initer:enditer; 
		figure,   
      Qpd=1/sigmap(kk).^2;
		Qp=Qpd*eye(np,np);      
      
      [V]=invforw(UH,Qp,CI); % Simplest inversion (Linear)
		
		V=V.';
		VD=duplic(V);
		vr=ifft(VD);
      
      [utr,JD(kk)]=backward(V,UH);
      
      JP(kk)=sum(dp*(sum(abs(V(1:60,:).^2))));
      J(kk)=JP(kk)+JD(kk);
      
      if Power==1 
         plotdata(ur,vr,utr,axis1)
      elseif Power==2
         plotdat2(ur,vr,utr,axis1)
      end   
      
      mytitle(sigmap(kk),pnorm);

end;
axis
plotnorm(sigmap,J,JP,JD);
mytitle2(pnorm);

