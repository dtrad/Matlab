function v = radonopid(d,t,h,p)
%RADONOPID Adjoint linear Radon operator
%
%   v = radonopid(d,t,h,p)
%
%   d  : data, size nt x nh
%   t  : time/tau axis
%   h  : offset axis
%   p  : Radon parameter axis
%
%   v  : Radon model, size nt x np
%
%   This is the adjoint of radonopd:
%
%       time = tau + p*h
%
%   Forward:
%       d(time,h) += v(tau,p)
%
%   Adjoint:
%       v(tau,p) += d(time,h)

nt = length(t);
nh = length(h);
np = length(p);
dt = t(2) - t(1);
v = zeros(nt,np,'like',d);
% if d is a vector reshape to d(nt,nh) and set flag=1
flag=0
if (isvector(d))
    d=reshape(d,nt,nh);
    flag=1
end

for ih = 1:nh
    hh = h(ih)^2;
    for ip = 1:np
        % Same shift as in the forward operator
        ishift = round(p(ip)*hh/dt);
        for itau = 1:nt
            it = itau + ishift;
            if it >= 1 && it <= nt
                v(itau,ip) = v(itau,ip) + d(it,ih);
            end
        end
    end
end
if (flag==1) 
    v = v(:);
end
end
