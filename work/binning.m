
function [x,index]=binning(xi,dx)
nx=length(xi);

j=1;
index(j)=1;
x(1)=round(xi(1)/dx)*dx;
for i=2:nx
    temp=round(xi(i)/dx)*dx;
    if (temp ~= x(j))
        j=j+1;
        x(j)=temp;
        index(j)=i;
    end
end
return;