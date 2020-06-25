function isitmp2(A,B)
% For A only use isitmp2(A)
% Plots zeros and poles to see whether A is 
% minimum phase (all zeros inside unit circle)
% maximum phase (all zeros outside unit circle)
% Note that convention here is similar to Electrical Ingeneering
% and different to Geophysics
% Daniel Trad -- UBC-- 6-08-98
if nargin < 2 B=[0];else B=B./B(1);end

zpplot(th2zp(poly2th(A./A(1),B )));
figure(gcf)