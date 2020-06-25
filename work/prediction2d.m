function [data3]=fxprediction
%example fxprediction interpolation for regular aliased data;
noise = 0.5;
over = 2;
gaps = 0;
load lineregnoalias
[nt nh]=size(data);

% Insert zero traces in original data;
zerotrace = zeros(nt,1);
data2=data(:,1);
for i=1:nh-1
    data2(:,2*i)=zerotrace;
    data2(:,2*i+1)=data(:,i+1);
end
data2(:,2*nh)=zerotrace;

% generate h2 regular from h 
dh = (h(end)-h(1))/(length(h)-1);
h2=h(1):dh/over:h(end);

datan=data+(rand(size(data))-0.5)*noise;

% make some traces equal to zero to see effect of gaps
if (gaps == 1) 
    datan(:,10)=0;
    datan(:,20)=0;
end

data3=filter(datan,t,h);

figure,wigb(datan);title('noisy')
figure,wigb(data3);title('clean');

return;


function [r]=autocor(d,l1,l2)
sum = 0;
for i1=1:n1
    for i2=1:n2
        if (((i1+l1)<= n1)&((i2+l2)<=n2))
            sum = sum + d(i1,i2)*d(i1+l1,i2+l2);
        end
    end
end


function [a]=prediction2d(d)
lag1=1:3;
lag2=1:3;

for k=1:M
   k1 = lag1(k);
   k2 = lag2(k);
   for j=1:M
       j1 = lag1(j);
       j2 = lag2(j);
       R(k,j)=autocor(d,k1-j1,k2-j2);
   end
   r(k)=autocor(d,k1,k2); 
end
a=R\r;
for j=1:M
   j1 = lag1(j);
   j2 = lag2(j);
   a(j1,j2)=a(j); 
end