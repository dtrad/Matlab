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