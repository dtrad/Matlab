
dz=10;
dx=10;
nz = 51;
nx = 51;
iz = 23;
vel=2500*ones(nz,nx);
for iii=1:nz
    if (iii<iz) 
        v=2000;
    else
        v=3000;
    end
    vel(iii,:)=v;
end
vel_initial=vel;
vel_initial(:,:)=2200;

for iii=1:nz
    for jjj=1:nx
        if ((iii-15)^2+(jjj-25)^2)<25
            vel(iii,jjj)=3000;
        end
        if ((iii-32)^2+(jjj-25)^2)<25
            vel(iii,jjj)=2200;
        end
    end
end

vel_true=vel;


