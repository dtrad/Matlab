
%Forwal modeling-elastic wave equation 

clear all
close all

dx=4;
dt=.002;
tmax= 0.76;
fdom=15;
dx0=50;
dz0=50;
dz=dx;
rho=1200;
qmax=1e-7;
n=2;
% save_snaps=1;

%--------Velocity model-------------------------------
   v=zeros(400,400);
   vP=zeros(size(v));
   vS=zeros(size(v));
   vP(:,:)=2500;
   vS(:,:)=1500;
    
%----Location of source and receivers-----------------  
[nz,nx] = size(v);
x=0:dx:(nx-1)*dx;
z=0:dx:(nz-1)*dx; 
xrec=0:1*dx:(size(v,2)-1)*dx;
zrec=xrec*0+3*dx;
zsor=200*dx;
xsor=(size(v,2)-1)*dx/2;

%--------delta t-----------------------------------------
dtstep = 0.0009;

%----Ricker Wavelet--------------------------------------
[wavelet,tw]=ricker(dt,fdom);
tw_2=min(tw):dtstep:max(tw); 
wavelet=interp1(tw,wavelet,tw_2);tw=tw_2;

%----Parameters Definition-------------------------------

nt=fix(tmax/dtstep)+1;

indr=sub2ind(size(v), fix(zrec/dz)+1, fix(xrec/dx)+1);
inds=sub2ind(size(v), fix(zsor/dz)+1, fix(xsor/dx)+1);


tic

t=0:dtstep:tmax;

%--PML Absorbing Boundary Condition----------------------

ind_mqx=0:dx0;ind_mqx=qmax*(ind_mqx/dx0).^n;
ind_mqz=0:dz0;ind_mqz=qmax*(ind_mqz/dz0).^n;

qpx=zeros(1,nx);qpz=zeros(1,nz);


qpx(1:length(ind_mqx))=fliplr(ind_mqx);qpx(end-length(ind_mqx)+1:end)=(ind_mqx);
qpz(1:length(ind_mqz))=fliplr(ind_mqz);qpz(end-length(ind_mqz)+1:end)=(ind_mqz);
qpz=qpz(:);
Qx=zeros(nz,nx);
Qz=zeros(nz,nx);

for k=1:nx
    Qz(:,k)=qpz; Qz(:,end-k+1)=qpz; 
end


for k=1:nz
    Qx(k,:)=qpx; Qx(end-k+1,:)=qpx;
end

if 1<0
Qz(1:dx0,:)=sqrt(Qx(1:dx0,:).*Qz(1:dx0,:));
end
QX=Qx;QZ=Qz;
pp=find(Qx==Qz);
Qx(pp)=QX(pp)+QZ(pp);
Qz(pp)=QX(pp)+QZ(pp);
Qx1=Qx;
% Qx=Qx*0;Qz=Qz*0;
clear QX QZ Qx1
figure;subplot(121);
imagesc(x,z,Qz)
axis equal
subplot(122);
imagesc(x,z,Qx)

%-----------Relaxed Modulus------------------------------
P1  = zeros(nz,nx,2);
P2= zeros(nz,nx,2);
P3= zeros(nz,nx,2);
Ux2  = zeros(nz,nx,2);
Uz2= zeros(nz,nx,2);
dataP1  = [];
dataP2 = [];
dataP3 = [];
dataUx2 = [];
dataUz2 = [];

Pdot=v*0;PPPdot=v*0;Uxdot=v*0;Uzdot=v*0;
PPdot=v*0;

zin = 4:(nz-3);     
xin = 4:(nx-3);     

x=0:dx:(nx-1)*dx;
z=0:dx:(nz-1)*dx; 
E_snap2=[];  % initialize 1D vector for ploting energy of P wavefield

kk=0;
kkk=0;   

M=vS./vP;
M11=rho*(vP.^2);
M33=rho*(vP.^2);
M44=rho*(vS.^2);
M13=-2.*M44+M33;   

%%
for it = 1:nt    
    
    if it<length(wavelet)   
        Uz2(inds)=Uz2(inds)+wavelet(it)*1;
    end
 
    
Ux2(zin,xin,2) = (2*rho-dtstep*(Qx(zin,xin)+Qz(zin,xin)))./(2*rho+dtstep*(Qx(zin,xin)+Qz(zin,xin))).*Ux2(zin,xin,1)...
               +2*dtstep./(2*rho+dtstep*(Qx(zin,xin)+Qz(zin,xin)))/24/dx.*((-P1(zin,xin+2-1,1)+27*P1(zin,xin+1-1,1)-27*P1(zin,xin-1,1)+P1(zin,xin-2,1))+(-P2(zin+2-1,xin,1)+27*P2(zin+1-1,xin,1)-27*P2(zin-1,xin,1)+P2(zin-2,xin,1)))...
               -(((Qx(zin,xin).*Qz(zin,xin)))).*Uxdot(zin,xin);

Uz2(zin,xin,2) = (2*rho-dtstep*(Qx(zin,xin)+Qz(zin,xin)))./(2*rho+dtstep*(Qx(zin,xin)+Qz(zin,xin))).*Uz2(zin,xin,1)...
                +2*dtstep./(2*rho+dtstep*(Qx(zin,xin)+Qz(zin,xin)))/24/dz.*((-P2(zin,xin+2-1,1)+27*P2(zin,xin+1-1,1)-27*P2(zin,xin-1,1)+P2(zin,xin-2,1))+(-P3(zin+2-1,xin,1)+27*P3(zin+1-1,xin,1)-27*P3(zin-1,xin,1)+P3(zin-2,xin,1)))...
                -(((Qx(zin,xin).*Qz(zin,xin)))).*Uzdot(zin,xin);
 
 P1(zin,xin,2) = (((2./M33(zin,xin))-dtstep*(Qx(zin,xin)+Qz(zin,xin)))./((2./M33(zin,xin))+dtstep*(Qx(zin,xin)+Qz(zin,xin)))).*P1(zin,xin,1)...
             +(((2*dtstep)./((2./M33(zin,xin))+dtstep*(Qx(zin,xin)+Qz(zin,xin)))))/24/dx.*(-(Ux2(zin,xin+2,2)+Qz(zin,xin+2).*Uxdot(zin,xin+2))+27*(Ux2(zin,xin+1,2)+Qz(zin,xin+1).*Uxdot(zin,xin+1))-27*(Ux2(zin,xin-1+1,2)+Qz(zin,xin-1+1).*Uxdot(zin,xin-1+1))+(Ux2(zin,xin-2+1,2)+Qz(zin,xin-2+1).*Uxdot(zin,xin-2+1)))...
             +(((2*dtstep)./((2./M33(zin,xin))+dtstep*(Qx(zin,xin)+Qz(zin,xin)))).*(-2.*(M(zin,xin).^2)+1))/24/dx.*(-(Uz2(zin+2,xin,2)+Qx(zin+2,xin).*Uzdot(zin+2,xin))+27*(Uz2(zin+1,xin,2)+Qx(zin+1,xin).*Uzdot(zin+1,xin))-27*(Uz2(zin-1+1,xin,2)+Qx(zin-1+1,xin).*Uzdot(zin-1+1,xin))+(Uz2(zin-2+1,xin,2)+Qx(zin-2+1,xin).*Uzdot(zin-2+1,xin)))...
             -(((Qx(zin,xin).*Qz(zin,xin)))).*Pdot(zin,xin);
        
 P3(zin,xin,2) = (((2./M33(zin,xin))-dtstep*(Qx(zin,xin)+Qz(zin,xin)))./((2./M33(zin,xin))+dtstep*(Qx(zin,xin)+Qz(zin,xin)))).*P3(zin,xin,1)...
             +(((2*dtstep)./((2./M33(zin,xin))+dtstep*(Qx(zin,xin)+Qz(zin,xin)))).*(-2.*(M(zin,xin).^2)+1))/24/dx.*(-(Ux2(zin,xin+2,2)+Qz(zin,xin+2).*Uxdot(zin,xin+2))+27*(Ux2(zin,xin+1,2)+Qz(zin,xin+1).*Uxdot(zin,xin+1))-27*(Ux2(zin,xin-1+1,2)+Qz(zin,xin-1+1).*Uxdot(zin,xin-1+1))+(Ux2(zin,xin-2+1,2)+Qz(zin,xin-2+1).*Uxdot(zin,xin-2+1)))...
             +(((2*dtstep)./((2./M33(zin,xin))+dtstep*(Qx(zin,xin)+Qz(zin,xin)))))/24/dx.*(-(Uz2(zin+2,xin,2)+Qx(zin+2,xin).*Uzdot(zin+2,xin))+27*(Uz2(zin+1,xin,2)+Qx(zin+1,xin).*Uzdot(zin+1,xin))-27*(Uz2(zin-1+1,xin,2)+Qx(zin-1+1,xin).*Uzdot(zin-1+1,xin))+(Uz2(zin-2+1,xin,2)+Qx(zin-2+1,xin).*Uzdot(zin-2+1,xin)))...
             -(((Qx(zin,xin).*Qz(zin,xin)))).*PPPdot(zin,xin);
              
        
P2(zin,xin,2) = (((2./M33(zin,xin))-dtstep*(Qx(zin,xin)+Qz(zin,xin)))./((2./M33(zin,xin))+dtstep*(Qx(zin,xin)+Qz(zin,xin)))).*P2(zin,xin,1)...
             +((((2*dtstep)./((2./M33(zin,xin))+dtstep*(Qx(zin,xin)+Qz(zin,xin)))).*(M(zin,xin).^2)))/24/dx.*((-(Uz2(zin,xin+2,2)+Qx(zin,xin+2).*Uzdot(zin,xin+2))+27*(Uz2(zin,xin+1,2)+Qx(zin,xin+1).*Uzdot(zin,xin+1))-27*(Uz2(zin,xin-1+1,2)+Qx(zin,xin-1+1).*Uzdot(zin,xin-1+1))+(Uz2(zin,xin-2+1,2)+Qx(zin,xin-2+1).*Uzdot(zin,xin-2+1)))+(-(Ux2(zin+2,xin,2)+Qz(zin+2,xin).*Uxdot(zin+2,xin))+27*(Ux2(zin+1,xin,2)+Qz(zin+1,xin).*Uxdot(zin+1,xin))-27*(Ux2(zin-1+1,xin,2)+Qz(zin-1+1,xin).*Uxdot(zin-1+1,xin))+(Ux2(zin-2+1,xin,2)+Qz(zin-2+1,xin).*Uxdot(zin-2+1,xin))))...
             -(((Qx(zin,xin).*Qz(zin,xin)))).*PPdot(zin,xin);

         
                                   
 Pdot=Pdot+P1(:,:,2);
 PPdot=PPdot+P2(:,:,2);
 PPPdot=PPPdot+P3(:,:,2);
 Uxdot=Uxdot+Ux2(:,:,2);
 Uzdot=Uzdot+Uz2(:,:,2); 
 

% ------------Iterations-------------------------------------
if it<nt   
 P1(:,:,1) = P1(:,:,2);
 P2(:,:,1) = P2(:,:,2);
 P3(:,:,1) = P3(:,:,2);
 Ux2(:,:,1) = Ux2(:,:,2);
 Uz2(:,:,1) = Uz2(:,:,2);

end

if mod(it*dtstep,dt)<dtstep
    kk=kk+1;
    dataUx2(kk,:) = Ux2(indr+nx*nz);
    dataUz2(kk,:) = Uz2(indr+nx*nz);
    dataP2(kk,:) = P2(indr+nx*nz);
    dataP3(kk,:) = P3(indr+nx*nz);
    dataP1(kk,:) = P1(indr+nx*nz);
end

if mod(it,10)==0


figure(100);

    imagesc(x,z,real(Uz2(:,:,2)))
       xlabel('Distance (m)','fontSize',25); ylabel('Depth (m)','fontSize',25);
       title(['Time = ',num2str(it*dtstep,'%10.3f')],'fontSize',25)
        colormap('seismic')
       axis equal 
       caxis([-1 1]*1e-2)




       kkk=kkk+1;
       E_snap2(kkk)=sum(sum(Uz2(:,:,2).^2));
       time_energy(kkk)=it*dtstep;



end
end

 figure(100);

    imagesc(x,z,real(Uz2(:,:,2)))
       xlabel('Distance (m)','fontSize',20); 
       ylabel('Depth (m)','fontSize',20);
       set(gca, 'FontSize', 15);
       caxis([-1 1]*1e-2)
       colormap('seismic')
       axis square
       print('WF', '-dpng', '-r400');
        
toc


