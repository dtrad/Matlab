function [f,g,Hv,imag3] = Inversion_vd_parameterization_encoding_pwy(model,D,frequency,nz,nx,dz,dx,xs,xr,sz,rz,w,w_f)
%%

        m=model(:,1).^2;
        density_inv=model(:,2);
        
    

%%
% get sampling operators
    Pr=zeros(size(xr,2),nx*nz);

    for ir=1:size(xr,2)
        Pr_temp=zeros(nz,nx);
        Pr_temp(rz,xr(1,ir))=1;
        Pr(ir,:)=reshape(Pr_temp,1,nx*nz);
    end
    clear Pr_temp;

    Ps=zeros(size(xs,2),nx*nz);
    for is=1:size(xs,2)
        Ps_temp=zeros(nz,nx);
        Ps_temp(sz,xs(1,is))=1;
        Ps(is,:)=reshape(Ps_temp,1,nx*nz);
    end
    clear Ps_temp;
    
% initialize misfit and gradient
f = 0;
g = zeros(size(model));
g_v=zeros(size(model(:,1)));
imag3=zeros(size(model(:,1)));
% loop over frequencies 
for k = 1:size(frequency,2)

          Lk = Helmoholtz2D_vd_parameterization2(frequency(k),m,density_inv,[dz dx],[nz,nx],model);   %% impedance matrix L use initial model
          Uk = Lk\(Ps');  %% wavefield calculate use initial model
          Sk = Pr*Uk;    %%  receiver data 
          % add source effect 
          f_in=find(w_f==frequency(k));     
          Sk=Sk*abs(w(f_in));           
          RES=(Pr'*(Sk-D(:,:,k)));          
          Vk = Lk'\RES; %% data residual as the source          
          clear Lk RES;
          %RES=(Sk-D(:,:,k));
          %figure;imagesc(abs(RES));colorbar;title('RES') ATA');
	      f = f + .5*norm(Sk - D(:,:,k),'fro')^2;
          clear Sk;
          dLk = getdL_vd_parameterization1(frequency(k),[nz,nx],m);   %derivative to L
         % g_v=g_v-2*sum(conj(Uk).*(dLk'*(Vk)),2).*(model(:,1).^3);   %% imaging condition
          g_v=g_v-2*sum(conj(Uk).*(dLk'*(Vk)),2).*(model(:,1).*model(:,2));   %% imaging condition
          %test simple imaging condition
          imag=(2*pi*frequency(k))^2*Uk.*conj(Vk);
          imag2=sum(imag,2);
end
imag3=imag3+abs(imag2);
g(:,1)=g(:,1)+real(g_v);
%returns handle to Hv product
Hv = @(x)GN_vec_lsm(x,model,frequency,nz,nx,dz,dx,xs,xr,sz,rz,g,D);
end


function y = GN_vec_lsm(x,model,frequency,nz,nx,dz,dx,xs,xr,sz,rz,g,D)


        m=model(:,1).^2;
        density_inv=model(:,2);
        density=1./model(:,2);
  
    
% get sampling operators
    Ps=zeros(size(xs,2),nx*nz);
    for is=1:size(xs,2)
        Ps_temp=zeros(nz,nx);
        Ps_temp(sz,xs(1,is))=1;
        Ps(is,:)=reshape(Ps_temp,1,nx*nz);
    end
    clear Ps_temp;
    Pr=zeros(size(xr,2),nx*nz);
    
    for ir=1:size(xr,2)
        Pr_temp=zeros(nz,nx);
        Pr_temp(rz,xr(1,ir))=1;
        Pr(ir,:)=reshape(Pr_temp,1,nx*nz);
    end
    clear Pr_temp;
y = 0;

% loop over frequencies 
for k = 1:size(frequency,2)

            Lk = Helmoholtz2D_vd_parameterization2(frequency(k),m,density_inv,[dz dx],[nz,nx],model);
            Uk = Lk\(Ps');
            dLk = getdL_variable_density_v_inv(frequency(k),[nz,nx],m,density);            
            dSk = (Pr*(Lk\(dLk*spdiags(x(1:length(x)),0,length(x),length(x))*Uk)));  % take vector x as the diagonal of the length(x)*length(x) martrix 
            g_bulk_bulk=sum(conj(Uk).*((dLk'*(Lk'\(Pr'*dSk)))),2);  %  sum over row
    
    y=y+real(g_bulk_bulk(:));
end
end