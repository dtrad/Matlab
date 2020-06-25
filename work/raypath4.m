clear all;
close all;
z=0:10:2000;
vo=1800;c=.6;nrays=20;
v=vo+c*z;
thetamin=5;thetamax=80;
deltheta=(thetamax-thetamin)/nrays;
xraypath=zeros(nrays,length(z));
for k=1:nrays
    theta=thetamin+(k-1)*deltheta;
    p=sin(pi*theta/180)/vo;
    pk(k)=p;
    x = (p*c)^-1 * ( sqrt(1-p^2*vo^2) - sqrt(1-p^2*(vo+c*z).^2) );
   
    ind=find(imag(x)~=0.0);
    if(~isempty(ind))
        x(ind)=nan*ones(size(ind));
    end
    ind=find(real(x)<0.);
    if(~isempty(ind))
        x(ind)=nan*ones(size(ind));
    end
    
    xraypath(k,:)=real(x);
    
end

%hold on;
plot(xraypath,z,'r',v,z,'g:');flipy;
%hold off

xsource=xraypath(:,end);

xlabel('x meters');ylabel('z meters')

save run1 pk xsource v;