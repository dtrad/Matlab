function [ Hout , gout , f ] = Hessian( frequency , fwave , vel , FDFD , D , R, nz , nx , PML_thick )
%%
%Define terms
f=0;
nzPML = nz + 2*PML_thick;
nxPML = nx + 2*PML_thick;
NN = nz*nx;
NPML = nzPML*nxPML;

for n=1:length(frequency)
    omega = 2*pi*frequency(n);
    
    [Gs,A] = FDFD(frequency(n),fwave(n),vel);

    if length(frequency) > 1
        res = R*Gs - D(:,:,n);
    else
        res = R*Gs - D;
    end
    
    Gr = A.'\R';
    
    %%
    gf = real(sum(-omega.^2*Gs.*(Gr*conj(res)),2));
    if n==1
        g = gf;
    else
        g = g + gf;
    end
    
    %%
    H_term = omega^4*real((Gr*Gr').*(Gs*Gs'));

    if n==1
        H = H_term;
    else
        H = H + H_term;
    end
    
    f = f+(norm(res,'fro')^2)/2;

end
%%
gout = zeros(NN,1);
Htemp = zeros(NPML,NN);
Hout = zeros(NN,NN);
for n=1:length(g)
    x_pos = floor(n/nzPML) + 1 - PML_thick;
    z_pos = n - (x_pos+PML_thick - 1)*nzPML - PML_thick;
    if z_pos < 1
        z_pos = 1;
    elseif z_pos > nz
        z_pos = nz;
    end
    if x_pos < 1
        x_pos = 1;
    elseif x_pos > nx
        x_pos = nx;
    end
    boundind = z_pos + nz*(x_pos-1);
    gout(boundind) = gout(boundind) + g(n);
    Htemp(:,boundind) = Htemp(:,boundind) + H(:,n);
end

for n=1:length(g)
    x_pos = floor(n/nzPML) + 1 - PML_thick;
    z_pos = n - (x_pos+PML_thick - 1)*nzPML - PML_thick;
    if z_pos < 1
        z_pos = 1;
    elseif z_pos > nz
        z_pos = nz;
    end
    if x_pos < 1
        x_pos = 1;
    elseif x_pos > nx
        x_pos = nx;
    end
    boundind = z_pos + nz*(x_pos-1);
    Hout(boundind,:) = Hout(boundind,:) + Htemp(n,:);
end
 


