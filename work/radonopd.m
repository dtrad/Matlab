function dr = radonopd(v,t,h,p)
% RADONOPD  Forward parabolic Radon transform (data <- model)
% dr = radonopd(v,t,h,p,adjoint)
% v(nt,np) -> dr(nt,nh). If v is a vector it returns dr as a vector.
% check dimensions of v are (nt,np)
% Ensure the input vector v has the correct dimensions

nt = length(t);
nh = length(h);
np = length(p);
flag=0;
% if v is a vector reshape to v(nt,np)
if isvector(v)
    v = reshape(v, nt, np);
    flag=1;
end

assert(size(v, 1) == nt && size(v, 2) == np, 'Input v must be of size (nt, np)');
dt = t(2) - t(1);

dr = zeros(nt,nh,'like',v);

for ih = 1:nh
    hh=h(ih)^2;
    for ip = 1:np
        % Time shift produced by p*h
        ishift = round(p(ip)*hh/dt);
        for itau = 1:nt
            it = itau + ishift;
            if it >= 1 && it <= nt
                dr(it,ih) = dr(it,ih) + v(itau,ip);
            end
        end
    end
end


if (flag==1) 
    dr=dr(:);
end
end


   