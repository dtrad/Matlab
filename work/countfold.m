function [fold,vnew,ix]=countfold(v,binsize)
if (nargin < 2) binsize=1;end;
[u,ix] = sort(round(v/binsize)*binsize);
nv = length(v);
% for speeding up, allocate vectors ahead and then ajust size
vnew = zeros(size(u));
fold = zeros(size(u));

i=2;
old = round(u(1));
inew = 1;
fold(inew) = 1;
vnew(inew)= u(1);
while (i<=nv)
    
    if (u(i) == old)
        fold(inew) = fold(inew) + 1;        
    else
        inew = inew + 1;
        vnew(inew)= u(i);
        fold(inew) = 1;
        old = u(i);
    end
    i = i + 1;
end

vnew = vnew(1:inew);
fold = fold(1:inew);

return 

