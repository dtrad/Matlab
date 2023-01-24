function dL = getdL_vd_parameterization1(f,n,m)
omega = 2*pi*f;

N     = prod(n);
w = [0 ones(1,n(1)-2) 0];
if n(2)>1
    w = w(:)*[0 ones(1,n(2)-2) 0];
end
w = w(:);

dL = omega^2*spdiags(w,0,N,N) ;