function [g_alpha,delta_x_alpha,image] = Optimization_LSRTM(Inversion,model_initial,model_true,nz,nx,cgmax)
%% Linear conjugate gradient LSRTM solving H*delta_m=-g  (Ax=b) 



%% initialization
fprintf('# ---------------------------------------------#\n');
fprintf('# Running LSRTM:\n');
converged = 1;
iter      = 0;


model_x=model_initial;

% alpha0    = 1;
cgiter    = 0;
fprintf('# entering inversion function \n');
% initial evaluation
[f,g,Hv,image] = Inversion(model_x);
fprintf('# leaving inversion function \n');
% main loop
% alpha0 = 1/norm(g(:));
if (converged==1)    
    g_temp=reshape(g,nx*nz,2);
    g_alpha(:)=g_temp(:,1);
    delta_x_alpha=g_alpha;
end
while ~converged

%        alpha0=1;
%        if iter==0
%             alpha0=1/(norm(g(:)));
%        end
       fprintf('# CG function \n');
       g_temp=reshape(g,nx*nz,2);
       g_alpha(:)=g_temp(:,1);
%      save g_alpha g_alpha;
       [delta_x_alpha,cgiter,FLAG]=CG(Hv,-g_alpha(:),zeros(size(g_alpha(:))),1e-3,cgmax);
%      save delta_x_alpha delta_x_alpha
       converged=1;
         fprintf('#Leaving CG function \n');   
end
    
    fprintf('## ---------------------------------------------##\n');
    fprintf('Inversion End.\n');   
end

