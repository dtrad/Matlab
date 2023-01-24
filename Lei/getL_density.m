function [L1st]= getL_density(h,n)
D1 = spdiags(ones(n(1),1)*[-1 1]/(2*h(1)),[-1:2:1],n(1),n(1));
D2 = spdiags(ones(n(2),1)*[-1 1]/(2*h(2)),[-1:2:1],n(2),n(2));
% D1 = spdiags(ones(n(1),1)*[-1 1]/(2*h(1)),[-1:2:1],n(1),n(1));
% D2 = spdiags(ones(n(2),1)*[-1 1]/(2*h(2)),[-1:2:1],n(2),n(2));
% D1(1,1:2)=[1 -1]./10;
% D1(end-1:end,end)=[-1 1]./10;
% D2(1,1:2)=[1 -1]./10;
% D2(end-1:end,end)=[-1 1]./10;
% save D1 D1;
% save D2 D2;

N=n(1)*n(2);

L_x=kron(D2,speye(n(1)));
L_z=kron(speye(n(2)),D1);

w = [0 ones(1,n(1)-2) 0];
if n(2)>1
    w = w(:)*[0 ones(1,n(2)-2) 0];
end
w = w(:);

W=spdiags([w w],[0:1],N,N);
% save L_x L_x;
% save L_z L_z;
% save W W;
% L_x = L_x.*W;
% L_z = L_z.*W;
%     
L1st=(L_x.*W+L_z.*W);
% L1st=L_x+L_z;



