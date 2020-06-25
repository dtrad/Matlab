function [z]=irregconv(y,x)
nx=length(x);
nz=2*nx-1;
z=zeros(1,nx);

for i=1:nz
    z(i)=0;
    for j=1:nx
        z(i)=z(i)+y(j)*func(x(j))
    end
end
    
        

return


function [y]=func(x)

y=exp(-x*x/10);

return