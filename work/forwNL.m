function [V,JPF,J,JD,JP]=ForwNL(UH,VLS,tol,NITER,tolmodel,pnorm,sigmap,Cn)
% Forward Transform (Non-linear): v=Lu or v=FWU.u
% u are the original data, 
% input 
%  	UH: t-x data
%		VLS: initial model
% 		tol: tolerance for small v(i)
%		NITER: number of iterations
%		tolmodel: minimum difference between iterations
%		pnorm: norm required 1 for L1, 2 for L2, 3 for Cauchy
%		sigmap: scale of the posterior distribution of parameters
%     CI: inverse of covariance of noise
% output 
% 		V:  w-p data
%     JPF: sum of v;
%		J: cost function
% 		JD: data misfit
%     JP: model norm 
% Daniel Trad-- 6-04-98
global w h0 nt Power WV WU0 W np alfa

h=h0;
WU=WU0;
nh= min(size(UH));

Cn=diag(diag(Cn));
CI=diag(1./diag(Cn));


ii=1:np;
JPF=zeros(nt/2,NITER);
J=zeros(nt/2,NITER);
JD=zeros(nt/2,NITER);
JP=zeros(nt/2,NITER);

for f=1:nt/2;
   %----------------------------------
   % Initial model
   vold0=zeros(size(VLS(f,ii))); 				% m_i-2
   vold1(ii)=zeros(size(VLS(f,ii))); 			% m_i-1
   v(ii)=vold1; 										% m_i
   Jold=1e20;
   %----------------------------------
   % Operators
   F=exp(i*w(f)*(alfa*(h.^Power)));
   FH=F';
   L=F*WU;
   LH=FH*WV;
   
   %----------------------------------
   for k=1:NITER;
   % Watch for changes on the model
  		toliter=sum(abs(v)).*tolmodel;
      v=v(:);vold0=vold0(:);vold1=vold1(:);
      if (sum(abs(v-vold1))<toliter)&(sum(abs(v-vold0))<toliter)&(k~=1)break,end
      vold0=vold1;vold1=v;
   % Huber condition to avoid NaN 
  		if pnorm==1 for iii=1:np; if abs(v(iii))<tol v(iii)=tol;end,end,end,
   % Qp matrix     
      if (pnorm==10)
         II=ones(size(ii)).';  
         Qpd(ii)=II./(sigmap^2)./(1+(abs(v(ii))./sigmap).^2)+1e-10;
      else
         Qpd(ii)=(1/2/sigmap^2).*(abs(v(ii))./sigmap).^(pnorm-2);
      end   
         Qp=diag(Qpd);
         % Two cases for inversion:
      if np < nh   
         v=inv(Qp+L*CI*LH)*L*CI*(UH(f,:)).';
      else    
      	QpI=diag(1./diag(Qp));
         v=QpI*L*inv(Cn+LH*QpI*L)*(UH(f,:)).';
      end
      ut=LH*v(:);
      
 		JD(f,k)=abs((W*(ut-UH(f,:).'))'*WU*(W*(ut-UH(f,:).')));     
     
      if pnorm==10
         JP(f,k)=diag(WV)'*log(1+(abs(v)./sigmap).^2);
      elseif pnorm==1
         JP(f,k)=diag(WV)'*((abs(v)./sigmap).^pnorm);
      end   
      
      J(f,k)=JP(f,k)+JD(f,k);
      
      if ((abs(J(f,k))*1.0./abs(Jold)) > 1 &(k~=1))
         v=vold1;
         display('cost function increases');
         f,k,
         break;
      end;
      Jold=J(f,k);
      
      if pnorm==11 maxvalue=max(abs(v));v=v/maxvalue;end;
      JPF(f,k)=sum(abs(v));
    

   end;
   if pnorm==11 v=v.*maxvalue;end;
	if(f==1) V=v;else V=[V,v];end,   
end,
