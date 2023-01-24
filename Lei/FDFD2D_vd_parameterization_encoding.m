function Data= FDFD2D_vd_parameterization_encoding(model,f,nz,nx,dz,dx,xs,xr,sz,rz,w,w_f)
%% model

m=model(:,1).^2;
density_inv=model(:,2);

%%
% define receiver vectors
Pr=zeros(size(xr,2),nx*nz);

% initialize Pr
    for ir=1:size(xr,2)
        Pr_temp=zeros(nz,nx);  % every receiver is nz*nx
        Pr_temp(rz,xr(1,ir))=1; % put 1 as the receiver exist
        Pr(ir,:)=reshape(Pr_temp,1,nx*nz);
    end
    
% initialize Ps
    Ps=zeros(size(xs,2),nx*nz);
    for is=1:size(xs,2)
        Ps_temp=zeros(nz,nx);
        Ps_temp(sz,xs(1,is))=1;
        Ps(is,:)=reshape(Ps_temp,1,nx*nz);
    end
    Data= zeros(size(xr,2),size(xs,2),size(f,2)); % Data is a 3D matrix, z axis indicates receivers, x axis indicates sources

% loop over frequencies  
for f_index = 1:size(f,2)    
	% get Helmholtz operator
    kkk=f_index;
	Lk = Helmoholtz2D_vd_parameterization2(f(kkk),m,density_inv,[dz dx],[nz,nx],model);      
%  save Lk;
    Uk = Lk\Ps';
  	Data(:,:,f_index) = Pr*Uk;  
    % add source effect
      f(f_index)
      f_in=find(w_f==f(f_index))     
      Data(:,:,f_index)=Data(:,:,f_index)*abs(w(f_in));
end
% save U U;
% Data=Data*w(10);

