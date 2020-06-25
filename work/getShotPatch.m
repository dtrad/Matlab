function [rpatchX,rpatchY,thisShotX,thisShotY,fold]=getShotPatch(shotID,shotx,shoty,recx,recy,oneID)

n = length(shotID);

% for speeding up, allocate vectors ahead and then ajust size
rpatch = zeros(size(shotID));
fold   = 0;

for i=1:n
    if (fold == 0)&(shotID(i) == oneID)
        thisShotX = shotx(i);
        thisShotY = shoty(i);
    end
    if (shotID(i) == oneID)
        fold = fold + 1;
        rpatchX(fold) = recx(i);
        rpatchY(fold) = recy(i); 
    end
end

rpatchX = rpatchX(1:fold);
rpatchY = rpatchY(1:fold);

return 

