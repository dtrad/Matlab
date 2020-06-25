function plotdat2(ur,vr,utr,axis1)
global dh nh h_near t dp np Power p alfa h

		subplot(221),
      wigb(real(ur),1,h,t);
      %plotrace(ur,dh,nh,h_near,t);
		title('original data in t-x domain'),xlabel('time'),ylabel('offset'),
      %axis(axis1);
      
      axis2=[axis1(1) axis1(2) alfa(1) alfa(np)];
      dp2=(alfa(np)-alfa(1))/np;
      
		subplot(222),
      temp=vr./max(max(abs(vr)));
      
      wigb(real(vr),1,alfa,t);
      %plotrace(temp,dp2,np,alfa(1),t);
      if Power==1
         title('data in tau-p domain'),xlabel('tau'),ylabel('p'),
      elseif Power==2
         title('data in parabolic tau-p domain'),xlabel('tau'),ylabel('p'),
		end
      %axis(axis2);
      
		subplot(223),
      wigb(real(utr),1,h,t);
      %plotrace(utr,dh,nh,h_near,t);
		title('recovered data in t-x domain'),xlabel('time'),ylabel('offset'),
      %axis(axis1);
      
		res=real(ur-utr);
      subplot(224),
      wigb(res,0.5,h,t);
      %plotrace(res,dh,nh,h_near,t);
		title('residuals between input and output with inverse problem')
		xlabel('time'),ylabel('offset'),
      %axis(axis1);    

