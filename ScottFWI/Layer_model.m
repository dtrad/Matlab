
dz=10;
dx=10;
nz = 40;
nx = 60;
vel=2500*ones(nz,nx);

vel_initial=vel;

for iii=1:nz
    for jjj=1:nx
        if iii > 10
            vel(iii,jjj)=2700;
        end
        if iii > 25
            vel(iii,jjj)=2900;
        end
    end
end

vel_true=vel;


