function dL = getdL_variable_density_v_inv(f,n,m,density)
omega = 2*pi*f;

v=1./sqrt(m(:));

N     = prod(n);
w = [0 ones(1,n(1)-2) 0];
if n(2)>1
    w = w(:)*[0 ones(1,n(2)-2) 0];
end
w = w(:);

dL = 2*omega^2*spdiags(w./(v(:).*density(:)),0,N,N);