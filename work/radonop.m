function dr = radonop(v,t,h,p)
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

% create output (same as before)
dr = zeros(nt, nh, 'like', v);

% parallel loop over receivers
parfor ih = 1:nh
    tempCol = zeros(nt,1,'like',v);   % local accumulator for this ih
    hh = h(ih)^2;
    for ip = 1:np
        ishift = round(p(ip)*hh/dt);
        for itau = 1:nt
            it = itau + ishift;
            if it >= 1 && it <= nt
                tempCol(it) = tempCol(it) + v(itau, ip);
            end
        end
    end
    dr(:, ih) = tempCol;  % single assignment back to sliced variable
end

if (flag==1) 
    dr=dr(:);
end
end


   