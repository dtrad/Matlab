nt = 501;
nh = 51;
np = 41;

dt = 0.004;

t = (0:nt-1)'*dt;
h = linspace(-1000,1000,nh)';
p = linspace(-0.001,0.001,np)';

% Random model and data
v = randn(nt,np);
d = randn(nt,nh);

% Forward operation
Lv = radonopd(v,t,h,p);

% Adjoint operation
Ltd = radonopid(d,t,h,p);

% Inner products
a = sum(Lv(:).*d(:));
b = sum(v(:).*Ltd(:));

fprintf('<Lv,d>   = %.15e\n',a);
fprintf('<v,L''d> = %.15e\n',b);
fprintf('relative error = %.15e\n',abs(a-b)/max(abs(a),abs(b)));