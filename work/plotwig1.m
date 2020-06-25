function plotwig1(ur,vr,utr,t,h0,h1,alfa,HH)
 
		      
		subplot(221),
      wigb(real(ur),1,h0,t);
      title('original data in t-x domain'),
      ylabel('time'),
      xlabel('offset'),
      vv=axis;

		subplot(222),
      vr=normalize(vr);
      wigb(real(vr),0.1,alfa,t);
      title('data in tau-p domain'),xlabel('tau'),ylabel('p'),
      
      
      subplot(223),
      wigb(real(utr),1,h1,t);
      title('recovered data in t-x domain'),
      ylabel('time'),
      xlabel('offset'),
                  
      utrtemp=shrinkt2(utr,HH);
      res=ur-utrtemp;
      subplot(224),
      amx=max(max(abs(ur)));
      wigb(real(res),1,h0,t,amx);
      title('residuals')
      ylabel('time'),
      xlabel('offset'),
       
		subplot(221);axis(vv);
      subplot(223);axis(vv);
      subplot(224);axis(vv);
