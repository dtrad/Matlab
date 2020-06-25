function plotwig2(ur,vr,utr,t,alfa,h0,h1,HH,residuals)

if isempty(residuals) residuals='n';end % Change to yes to produce residual plot
 
		subplot(221),
      wigb(real(ur),1,h0,t);
      title('(a) Original data in the t-x domain'),
      ylabel('time'),xlabel('offset'),
           
		subplot(222),
      wigb(real(vr),1,alfa,t);
      title('(b) Data in the Radon domain'),
      ylabel('tau'),xlabel('q'),
      
		subplot(223),
      wigb(real(utr),1,h1,t);
      title('(c) Recovered data in the t-x domain'),
      ylabel('time'),
      xlabel('offset'),
      vv=axis;
      
      if residuals=='y'
      utrtemp=shrinkt2(utr,HH);
		res=real(ur-utrtemp);
      amx=max(max(abs(ur)));
      subplot(224),
      wigb(res,1,h0,t,amx);
      title('(d) Residuals')
      ylabel('time'),
      xlabel('offset'),
      subplot(224);axis(vv);
      end
   
		subplot(221);axis(vv);
      subplot(223);axis(vv);
         