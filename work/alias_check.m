function [message]=alias_check(p,h,fmax)
% Check for possible alias
% [message]=aliastaup(p,h,fmax)
% Daniel Trad - UBC- 27-08-98
if nargin < 3 fmax=1/2/0.004;end
message=[];
dp=abs(p(1)-p(2))
deltah=abs(h(length(h))-h(1));
dpmax=abs(1/fmax/2/deltah)
if (dp>dpmax) message='possible alias';display(message);end



