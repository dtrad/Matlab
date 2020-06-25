function [shiftOffset,shiftAzimuth] = offsetShift(ix,iy,ia,ih,deltaH,deltaA)
% calculate an offset shift that can be dependent on any of cmpxy, azimuth
% or offset
ix,iy,ia,ia

% offset shift depending on azimuth
if (mod(ia,2)) 
    shiftOffset = deltaH;
else
    shiftOffset = 0;
end

shiftAzimuth = 0;
return;