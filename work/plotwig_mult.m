function plotwig_mult(ur,vrorig,utr,vr,t,alfa,h0,h1)
dt=0.004;
nh0=length(h0);
nt=length(t);
nh1=length(h1);
nalfa=length(alfa);
tau=t;
%tau=t.^2;
figure,
subplot(221);
      wigb(real(ur(1:nt,1:nh0)),1,h0,t);
      title('(a) Original data in the t-x domain'),
      ylabel('time'),xlabel('offset'),
           
subplot(222);
      wigb(real(vrorig(1:nt,1:nalfa)),1,alfa,tau);
      title('(b) Data in the Radon domain'),
      ylabel('tau'),xlabel('q'),
      
subplot(223);
      wigb(real(utr(1:nt,1:nh1)),1,h1,t);
      title('(c) Recovered data in the t-x domain'),
      ylabel('time'),
      xlabel('offset'),
      vv=axis;
     
   
subplot(224);
      wigb(real(vr(1:nt,1:nalfa)),1,alfa,tau);
      title('(d) Data in the Radon domain'),
      ylabel('tau'),xlabel('q'),
