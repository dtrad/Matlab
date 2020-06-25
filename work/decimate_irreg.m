function [dd,xx,yy,count,ds]=decimate_irreg(d,x,y,nx,ny,threshold)
count=0;
ds = zeros(size(d));
for ix=1:nx 
    for iy=1:ny
        if (rand() > threshold)
            count=count+1;
            dd(count)=d(ix,iy);
            xx(count)=x(ix);
            yy(count)=y(iy);
            ds(ix,iy)=1;
        end
    end
end
xx=xx(:);
yy=yy(:);
dd=dd(:);

return;