function v = radonopi(d,t,h,p)
%RADONOPID Adjoint linear Radon operator
nt = length(t);
nh = length(h);
np = length(p);
dt = t(2) - t(1);

flag = 0;
if isvector(d)
    d = reshape(d, nt, nh);
    flag = 1;
end

% allocate per-ih contributions (sliced on 3rd dim)
temp = zeros(nt, np, nh, 'like', d);

parfor ih = 1:nh
    hh = h(ih)^2;
    tmp = zeros(nt, np, 'like', d);  % local nt x np accumulator
    for ip = 1:np
        ishift = round(p(ip)*hh/dt);
        for itau = 1:nt
            it = itau + ishift;
            if it >= 1 && it <= nt
                tmp(itau, ip) = tmp(itau, ip) + d(it, ih);
            end
        end
    end
    temp(:, :, ih) = tmp;  % sliced assignment (safe)
end

% combine contributions from all ih
v = sum(temp, 3);

if flag
    v = v(:);
end
end
