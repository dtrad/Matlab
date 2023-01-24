function [rf,rg,rH] = Diff0_Reg(model,regfac,nz,nx)
  %%  
    
  rg=zeros(size(model));
  
  b1=zeros(1,2*nz*nx);
  b1(1:(nz*(nx-1)))=-1;
  b1(1+nz*nx:nz*nx+(nz*(nx-1)))=-1;
  
  b2=zeros(1,2*nz*nx);
  b2(nz+1:(nz*(nx)))=1;
  b2(1+nz+nz*nx:end)=1;
  
  b3=ones(1,2*nz*nx);
  b3(1:nz:nz*nx)=0;
  b3(1+nz*nx:nz:end)=0;
  
  b4=-1*ones(1,2*nz*nx);
  b4(nz:nz:nz*nx)=0;
  b4(nz+nz*nx:nz:end)=0;
  
  B1=[b1',b2'];
  d1=[-nz,0];
  
  derivteststorex=spdiags(B1,d1,length(model),length(model));
  
  B2=[b4',b3'];
  d2=[-1,0];
  
  derivteststorey=spdiags(B2,d2,length(model),length(model));
  
  rH=regfac*(derivteststorey'*derivteststorey + derivteststorex'*derivteststorex);
  rf=0;
        
    
end
