function [x]=dummyfunction(t,d)
% function [x]=dummyfunction(t,d) 
% Use this function together with quad to get the area of a vector d
% quad('dummyfunction',1,length(d),[],1,d) 
% Daniel Trad
%d=[
%    -0.0010 -0.0051   -0.0210   -0.0688   -0.1749   -0.3337   -0.4449   -0.3194 ...
%    0.1418    0.7272    1.0000    0.7272    0.1418   -0.3194   -0.4449...
%   -0.3337   -0.1749   -0.0688   -0.0210   -0.0051   -0.0010 ];


x=interp1(1:length(d),d,t,'spline');
%x=abs(x);
