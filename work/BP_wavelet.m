function [w,nw]=BP_wavelet(dt,f)
nf = 256
if1 = round(f(1)*nf*dt+1.1); 
if2 = round(f(2)*nf*dt+1.1); 
if3 = round(f(3)*nf*dt+1.1); 
if4 = round(f(4)*nf*dt+1.1);

for i=1:if1-1
  w(i) = 0.;
end

for i=if1:if2
  w(i) = (i - if1)/(if2-if1);
end

for i=if2+1:if3-1
  w(i) = 1; 
end

for i=if3:if4
  w(i) = 1.-(i-if3)/(if4-if3);
end

for i=if4+1:nf/2+1
  w(i) = 0.;
end

for i = nf/2+2:nf
  w(i) = w(nf+2-i);
end


for i = 1:nf
  aux(i) = w(i)+j*0.;
end

aux2=ifft(aux);
aux=aux2;

time = 2./(f(3)-f(2));
L = floor(time/dt); 


for i=1:L+1
  w(i) = real(aux(L-i+2));
end

for i=L+2:2*L+1
  w(i) = real(aux(i-L));
end

nw = 2*L+1;
pi2 = 8.*atan(1.);
for i = 1:nw
  w(i) = w(i)* (0.54-0.46*cos(pi2*(i-1)/(nw-1)));
end


