function plotdat3(ur,vr,utr,axis1)
global dh nh h_near t dp np Power h p
		urtemp=ur./max(max(abs(ur)));
      
		subplot(221),
      wigb(real(ur),1,h,t);
      %plotrace(urtemp,dh,nh,h_near,t);
		title('original data in t-x domain'),xlabel('time'),ylabel('offset'),
      %axis(axis1);
     
		subplot(222),
      vrtemp=vr./max(max(abs(vr)));
      wigb(real(vr),1,p,t);
      %plotrace(vrtemp,dp,np,0,t);
      if Power==1
         title('data in tau-p domain'),xlabel('tau'),ylabel('p'),
      elseif Power==2
         title('data in parabolic tau-p domain'),xlabel('tau'),ylabel('p'),
		end
      %axis([axis1(1) axis1(2) 0 np*dp]);
      
      utrtemp=utr./max(max(abs(utr)));
      subplot(223),
      wigb(real(utr),1,h,t);
		%plotrace(utrtemp,dh,nh,h_near,t);
		title('recovered data in t-x domain'),xlabel('time'),ylabel('offset'),
      %axis(axis1);
      
		res=urtemp-utrtemp;
		subplot(224),
      wigb(real(res),1,h,t);
      %plotrace(res,dh,nh,h_near,t);
		title('residuals between input and output with inverse problem')
		xlabel('time'),ylabel('offset'),
      %axis(axis1);    

