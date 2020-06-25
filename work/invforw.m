function [V]=invforw(UH,Qp,Cn)
% Forward Transform through the inverse: 
% v=inv(Qp+L.CI.L*).L.CI.u
% u are the original data, 
% input 
%		UH w-x data
% 		Qp regularization term
%     CI inverse of noise covariance matrix
% output 
%		V  w-p data 

% Daniel Trad-- 6-04-98

global w h0 nt np Power WV WU0 alfa

h=h0;
WU=WU0;

nh=min(size(UH));
Cn=diag(diag(Cn));
CI=diag(1./diag(Cn));

for f=1:nt/2;
   % Operators
   F=exp(i*w(f)*(alfa*(h.^Power)));
   FH=F';
   L=F*WU;
   LH=FH*WV;
   if np < nh   
         v=inv(Qp+L*CI*LH)*L*CI*(UH(f,:)).';
   else    
      	QpI=diag(1./diag(Qp));
         v=QpI*L*inv(Cn+LH*QpI*L)*(UH(f,:)).';
   end
   if(f==1) V=v;else V=[V,v];end,
end;
