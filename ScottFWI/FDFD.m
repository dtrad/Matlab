function [U,A] = FDFD(vel,vel0,frequency,S,fwave,PML_thick,R_theor,nz,nx,dz,dx)
%% This is a quick attempt
nzPML = nz + 2*PML_thick;
nxPML = nx + 2*PML_thick;
if length(frequency)>1
    U = zeros(nzPML*nxPML,size(S,2),length(frequency));
else
    U = zeros(nzPML*nxPML,size(S,2));
end

for n=1:length(frequency)
    omega = 2*pi*frequency(n);
    [ A ] = Make_Helmholtz(vel,vel0,omega,PML_thick,nz,nx,dz,dx,R_theor);
    u = A\(S*fwave(n));
    if length(frequency)>1
        U(:,:,n) = u;
    else
        U=u;
    end    
end

