function dr = radonop(v,t,h,p)

nt = length(t);
nh = length(h);
np = length(p);

dt = t(2) - t(1);

dr = zeros(nt,nh,'like',v);

for ih = 1:nh
    for ip = 1:np
        % Time shift produced by p*h
        ishift = round(p(ip)*h(ih)/dt);
        for itau = 1:nt
            it = itau + ishift;
            if it >= 1 && it <= nt
                dr(it,ih) = dr(it,ih) + v(itau,ip);
            end
        end
    end
end
end