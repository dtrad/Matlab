% Estimation of the transfer function of the earth with the S-stransform

load /home/dtrad/MT/datos.dat

dt=2;
nfreq=68;

t=0:max(size(datos))-1;
t=t*dt;

ex=datos(:,1);
ey=datos(:,2);
hx=datos(:,3);
hy=datos(:,4);


[Sex,f]=stransform2(ex,t,nfreq,fmin,fmax);
[Sey,f]=stransform2(ey,t,nfreq,fmin,fmax);
[Shx,f]=stransform2(hx,t,nfreq,fmin,fmax);
[Shy,f]=stransform2(hy,t,nfreq,fmin,fmax);


%Autopowers
cexex=conj(Sex).*Sex;
ceyey=conj(Sey).*Sey;
chxhx=conj(Shx).*Shx;
chyhy=conj(Shy).*Shy;

%Crosspowers
cexey=conj(Sex).*Sey;
cexhx=conj(Sex).*Shx;
cexhy=conj(Sex).*Shy;

chxex=conj(Shx).*Sex;
chxey=conj(Shx).*Sey;
chxhy=conj(Shx).*Shy;

ceyhy=conj(Sey).*Shy;

%Transfer function. One per time, at every frequency.
for j=1:nt
  for k=1:nf
    
    % check
    z1=inv([Shxhx(k,j) Shxhy(k,j);Shyhx(k,j) Shyhy(k,j)])*Shxex(k,j) ...
       Shyex(k,j)]);
    
    z2=inv([Shxhx(k,j) Shxhy(k,j);Shxhy(k,j) Shyhy(k,j)])*Sexhx(k,j) ...
       Sexhy(k,j)]);
    
    %check
    zz=[z1(:).'; z2(:).'];
    zzz(j,k,:,:)=zz;
  end
end

figure(1);
subplot(221);imagesc(abs(zzz(:,:,1,1)));  %zxx
subplot(222);imagesc(abs(zzz(:,:,1,2)));  %zxy
subplot(223);imagesc(abs(zzz(:,:,2,1)));  %zyx
subplot(224);imagesc(abs(zzz(:,:,2,2)));  %zyy







