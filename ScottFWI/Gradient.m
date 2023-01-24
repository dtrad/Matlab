function [f,gout,Gs_out] = Gradient( frequency , fwave , vel , FDFD , D , R , nz , nx , PML_thick , bound )
%%
if nargin < 10
    bound = 1;
end

f=0;
nzPML = nz + 2*PML_thick;
nxPML = nx + 2*PML_thick;
NN = nz*nx;

for n=1:length(frequency)
    omega = 2*pi*frequency(n);
    
    [Gs,A] = FDFD(frequency(n),fwave(n),vel);
    
    if size(D,3)~=1
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
    
    f = f+(norm(res,'fro')^2)/2;
    
    if n==1
       Gs_out = zeros(size(Gs,1),size(Gs,2),length(frequency)); 
    end
    
    Gs_out(:,:,n) = Gs;  
end
%%
gout = zeros(NN,1);

if bound == 1
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
    end
else
    g = reshape(g,nzPML,nxPML);
    gout = g(PML_thick+1:PML_thick+nz, PML_thick+1:PML_thick+nx);
    gout = gout(:);
end
        


