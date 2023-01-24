function dL = getdL_variable_density_den_inv(f,n,m,density)
omega = 2*pi*f;
v=1./sqrt(m(:));
N     = prod(n);
w = [0 ones(1,n(1)-2) 0];
if n(2)>1
    w = w(:)*[0 ones(1,n(2)-2) 0];
end
w = w(:);

dL = omega^2*spdiags(w.*m,0,N,N);
% dL = omega^2*spdiags(w./(v(:).^2),0,N,N);