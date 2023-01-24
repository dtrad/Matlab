function H = Helmoholtz2D_vd_parameterization2(f,m,density_inv,h,n,model)
omega = 2*pi*f;

% m=1./(bulk./density);
density=1./density_inv;
bulk=1./(m.*density_inv);

N     = n(1)*n(2);

% variable density
C3_1=zeros(N,1);

C3_1(1:N-1)=density(2:N);
C3_2=density;
C3=0.5*(1./C3_1+1./C3_2)./h(1).^2;

C2_1=zeros(N,1);
C2_1(2:N)=density(1:N-1);
C2_2=density;
C2=0.5*(1./C2_1+1./C2_2)./h(1).^2;

D1_variable_density=spdiags([C3 -C3-C2 C2],[-1:1],N,N); 

D1_mask=spdiags(ones(n(1),1)*[1 1 1],[-1:1],n(1),n(1)); 

D11_mask=kron(speye(n(2)),D1_mask);
D1_variable_density2=D11_mask.*D1_variable_density;

D1_mask3=spdiags(ones(n(1),1)*[0 0 0],[-1:1],n(1),n(1)); 
D1_mask3(1,1:2) = [1.2 -1.1]/h(1); 
D1_mask3(end,end-1:end) = [-1.1 1.2]/h(1);
D11_mask3=kron(speye(n(2)),D1_mask3);

D1_variable_density3=D1_variable_density2+D11_mask3;
D1_variable_density3(1,1)=0.1;
D1_variable_density3(end,end)=0.1;
% 
D1_variable_density=D1_variable_density3;
% z direction

C5_1=zeros(N,1);

C5_1(1:end-(n(1)))=density(n(1)+1:end);
C5_2=density;
C5=0.5*(1./C5_1+1./C5_2)./h(2).^2;

C4_1=zeros(N,1);
C4_2=density;
C4_1(n(1)+1:end)=density(1:end-(n(1)));

C4=0.5*(1./C4_1+1./C4_2)./h(2).^2;

C6=-C4-C5;

C5(end-2*n(1)+1:end)=-0.1;
C4(1:2*n(1))=-0.1;
C6(1:n(1))=0.1;
C6(end-n(1)+1:end)=0.1;
D2_variable_density=spdiags([C5 C6 C4],[-n(1):n(1):n(1)],N,N); 

%%
S = D1_variable_density+ D2_variable_density;

w = [0 ones(1,n(1)-2) 0];
w = w(:)*[0 ones(1,n(2)-2) 0];
w = w(:);

M = omega^2*spdiags(w./bulk,0,N,N) + 1i*omega*spdiags((1-w).*sqrt(m),0,N,N);

H = M + S;






