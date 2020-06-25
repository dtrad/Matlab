% Example of Null space and precondtioning
% Steepest Descent example for null space and preconditioning
% paper from Dave Nichols, SEP report N82, pag 186
% Daniel Trad UBC

B=[1 1];
d=[1];
[u,s,v]=svd(B'*B)
m=[0;0]
r=d-B*m
g=B'*r
alfa=fmin('residual',0,10,[],d,B,g)
m=m+alfa*g
d-B*m
pause
% Right Preconditioner

m=[0;0];
Mr=[10 0;0 1];
[u,s,v]=svd(Mr'*B'*B*Mr)
pause
mr=Mr\m;
g=Mr'*B'*d;
alfa=fmin('residual',0,10,[],d,B*Mr,g)
mr=mr+alfa*g
m=Mr*mr



