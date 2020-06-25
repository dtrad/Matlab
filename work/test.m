function test(n,factor)
v=[0,1000,2000,3000,4000,5000];
t=0:n-1;
v=t*1000;
subplot(211),plot(t,v,'o');
for i=1:n
  if (t(i)>2) 
    v(i)=v(i)*(1-factor*sin(i*pi/(2*n)));
  end;
end;

subplot(212),plot(t,v,'o',interp(t,10,2),interp(v,10,2));figure(gcf);

