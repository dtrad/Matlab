x=-50:.1:50;
s = 3;
k0 = ;
sumg = zeros(size(x));
sumg = [sumg sumg];
n = length(sumg);
k = k0/sqrt(pi*s);
for i=1:20
  xc=i*k0;
  g=k*exp(-(x-xc).^2/s);
  if (i==10) g0=g;end
  figure
  plot(g);
  zpad = n - length(g);
  sumg = sumg + [g zeros(1,zpad)];
end
figure;plot(sumg)
