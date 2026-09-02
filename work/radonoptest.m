clear
close all
nt = 1001;
t = linspace(0,1,nt)';

h = linspace(-1000,1000,101)';
p = linspace(-0.0005,0.0005,101)';

v = zeros(length(t),length(p));

% Spike at tau = 0.5 s and p = 0.0003 s/m
[~,itau] = min(abs(t-0.5));
[~,ip]   = min(abs(p-0.0003));

v(itau,ip) = 1;

dr = radonop(v,t,h,p);

imagesc(h,t,dr);
axis xy;
xlabel('Offset h');
ylabel('Time t');