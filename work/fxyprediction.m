function [data3,a]=fxyprediction
%example fxprediction interpolation for regular aliased data;
noise = 0.5;
over = 2;
gaps = 0;
M=3;
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
[a]=fxyfilter(datan,M);
%data3=filter(datan,t,h);
data3 = conv2(datan,a);
figure,wigb(datan);title('noisy')
figure,wigb(data3);title('clean');

return;


function [sum]=autocor2d(d,l1,l2)
sum = 0;
[n1 n2]=size(d);
for i1=1:n1
    for i2=1:n2
        il1=i1+l1;
        il2=i2+l2;
        if ( (il1<= n1) & (il2<=n2)& (il1 > 0) & (il2 > 0) )
            sum = sum + d(i1,i2)*d(i1+l1,i2+l2);
        end
    end
end
return

function [a2]=fxyfilter(d,M)
k=1;
for i1 =1:M
    for i2 =1:M
        if ((i1-1)~=0) && ((i2-1)~=0) 
            lag(k)   = i1 - 1;
            lag(k+1) = i2 - 1;
        end
        
        k = k + 2;
    end
end

lag
neq = length(lag)/2
for k=1:neq
   k1 = lag(2*(k-1)+1);
   k2 = lag(2*(k-1)+2);
   for j=1:1:neq
       if ((j1)==0) && ((j2)==0) w(k) = 0; 
       else w(k) = 1;
       end
       j1 = lag(2*(j-1)+1);
       j2 = lag(2*(j-1)+2);
       R(k,j)=w(k)*autocor2d(d,k1-j1,k2-j2);
   end
   [k1 k2]
   r(k)=autocor2d(d,k1,k2); 
end

a=R\r(:);

for j1=1:M
    for j2=1:M
        a2(j1,j2)=a(2*(j1-1)+j2);
    end
end
