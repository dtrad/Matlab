%create offset to angle matrix
z=0:1:2500;
h=0:1:4000;
nz=length(z);
nh=length(h);
map=zeros(nz,nh);
rd=400;
for i=1:nz
        for j=1:nh
            map(i,j)=offDepth2angle(z(i),rd,h(j));
        end
end
image(map)
colorbar
xlabel('offset');
ylabel('depth');
title('angle for OBN with receivers at 400 depth');
%axisx(h);
%axisy(z);