clear all;
close all;
z=0:10:2000;
vo=1800;c=.6;nrays=20;
thetamin=5;thetamax=80;
deltheta=(thetamax-thetamin)/nrays;
xraypath=zeros(nrays,length(z));
for k=1:nrays
    theta=thetamin+(k-1)*deltheta;
    p=sin(pi*theta/180)/vo;
    pk(k)=p;
    x = (p*c)^-1 * ( sqrt(1-p^2*vo^2) - sqrt(1-p^2*(vo+c*z).^2) );
    x2 = (p*c)^-1 * ( sqrt(1-p^2*vo^2) + sqrt(1-p^2*(vo+c*z).^2) );
   %t=c^-1*log( ((v0+z*c)/v0) .*  ((1+sqrt(1-p^2*v0^2))./(1+sqrt(1-p^2(v0+c.*z^2)))) ); 
   
    ind=find(imag(x)~=0.0);
    if(~isempty(ind))
        x(ind)=nan*ones(size(ind));
    end
    ind=find(real(x)<0.);
    if(~isempty(ind))
        x(ind)=nan*ones(size(ind));
    end
   
   
    ind=find(imag(x2)~=0.0);
    if(~isempty(ind))
        x2(ind)=nan*ones(size(ind));
    end
    ind=find(real(x2)<0.);
    if(~isempty(ind))
        x2(ind)=nan*ones(size(ind));
    end
    
    xraypath(k,:)=real(x);
    xraypath2(k,:)=real(x2);
    
end

hold on;
plot(xraypath,z);
%plot(xraypath2,z);
flipy;
hold off

xsource=xraypath(:,end);

xlabel('x meters');ylabel('z meters')

save run1 pk xsource;