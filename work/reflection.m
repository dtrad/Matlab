load run1;
%sort out only the rays which have reached the reflector 
ind=find(real(xsource>0));
sources=zeros(size(ind));
sources(1:max(ind)) = xsource(1:max(ind));

%trace the reflected rays back to the surface
z=2000:-10:0;
vo2=max(v);c2=-.6;nrays=length(sources);
v2=1*(vo2+c2*z);
xraypath2=zeros(nrays,length(z));
for k=1:nrays
    p=pk(k);
    x = (p*c2)^-1 * ( sqrt(1-p^2*vo2^2) - sqrt(1-p^2*(vo2+c2*z).^2) );
   %t=c^-1*log( ((v0+z*c)/v0) .*  ((1+sqrt(1-p^2*v0^2))./(1+sqrt(1-p^2(v0+c.*z^2)))) ); 
   
    ind=find(imag(x)~=0.0);
    if(~isempty(ind))
        x(ind)=nan*ones(size(ind));
    end
    ind=find(real(x)<0.);
    if(~isempty(ind))
        x(ind)=nan*ones(size(ind));
    end
    
    xraypath2(k,:)=real(x);
    
end

%shift the reflected rays to the appropriate sources
rays=zeros(nrays,length(z));
for i = 1:nrays,
    rays(i,:) = xraypath2(i,:)+sources(i);
end

hold on;
plot(rays,z,'r',v2,z,'g:');
hold off

xlabel('x meters');ylabel('z meters')