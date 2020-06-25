function [shiftOffset,shiftAzimuth] = getShift(icmp,ia,ih,deltaH,deltaA)
% calculate an offset shift that can be dependent on any of cmpxy, azimuth
% or offset


% offset shift depending on azimuth
if (mod(ia+icmp,2)) 
    shiftOffset = deltaH;
else
    shiftOffset = 0;
end

shiftAzimuth = 0;
return;