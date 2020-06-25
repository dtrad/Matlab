function [data3,datan,a,lag1,lag2]=fxyprediction(M)
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
[a,lag1,lag2]=fxyfilter(datan,M);
%data3=filter(datan,t,h);
%data3 = conv2(datan,a);
data3=applyfilter(datan,a,lag1,lag2);
figure,wigb(datan,1);title('noisy')
figure,wigb(data3,1);title('clean');
figure,wigb(datan-data3,1);title('noise');

figure
t=1:551;
hold
const =0;
for h1=1:nh;
    const = const + 1;
    plot(t,datan(:,h1)+const,t,data3(:,h1)+const,t,datan(:,h1)-data3(:,h1)+const)
end

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

function [p]=applyfilter(d,a,lag1,lag2);
M=length(a);
p=zeros(size(d));
[n m]=size(d);
for i1=1:n
    for i2=1:m
        for j =1:M
            j1 = lag1(j);
            j2 = lag2(j);           
            if (((i1-j1) >0)&&((i2-j2)>0))
                if (((i1-j1) <= n)&&((i2-j2) <= m))
                    p(i1,i2)=p(i1,i2)+a(j)*d(i1-j1,i2-j2);
                end
            end
        end
        
    end
end
return


function [a,lag1,lag2]=fxyfilter(d,M)
k=1
for i=-M:M
    for j=-M:M
        if (i~=0)||(j~=0)
            lag1(k)=i;
            lag2(k)=j;
            k=k+1;
        end
    end
end

%lag1 = [ 1 0 1 -1  0 -1  1 -1 ];
%lag2 = [ 0 1 1 -1 -1  0 -1  1] ;

[lag1(:) lag2(:)]

M = length(lag1);

for k=1:M
   k1 = lag1(k);
   k2 = lag2(k);
   for j=1:1:M
       j1 = lag1(j);
       j2 = lag2(j);
       R(k,j)=autocor2d(d,k1-j1,k2-j2);
   end
   
   r(k)=autocor2d(d,k1,k2); 
end

a=R\r(:);

return
