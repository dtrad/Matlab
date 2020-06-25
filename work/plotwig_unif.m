function plotwig2(ur,vr,utr,t,alfa,h0,h1)
dt=0.004;
nh0=length(h0);
nt=length(t);
nh1=length(h1);
nalfa=length(alfa);
tau=t;
%tau=t.^2;

		figure,
      wigb(real(ur(1:nt,1:nh0)),1,h0,t);
      title('(a) Original data in the t-x domain'),
      ylabel('time'),xlabel('offset'),
           
		figure,
      wigb(real(vr(1:nt,1:nalfa)),1,alfa,tau);
      title('(b) Data in the Radon domain'),
      ylabel('tau'),xlabel('q'),
      
		figure,
      wigb(real(utr(1:nt,1:nh1)),1,h1,t);
      title('(c) Recovered data in the t-x domain'),
      ylabel('time'),
      xlabel('offset'),
      vv=axis;
     
   
   	figure,
      subplot(211);plot(h0,'+');
      title('(a) original offset');
      xlabel('# of trace');
      ylabel('offset (m)');
      subplot(212);plot(h1,'+');
      title('(b) final offset')
      xlabel('# of trace');
      ylabel('offset (m)');