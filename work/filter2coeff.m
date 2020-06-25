function [c]=filter2coeff(a,n)
% [c]=filter2coeff(a,n)
c(1)=-1;R(1)=1;a(1)=1;V(1)=1;
i=1:n+1;c(i)=a(i);
for k=1:n
   j=n-k+2;   
   temp=(1-c(j)*conj(c(j)));
   if temp<(1e-10) temp=1e-10;end
   al=1/temp;
   be=c(j)*al;
   jh=(j+1)/2;
   for i=1:jh
      top=al*c(i)-be*conj(c(j-i+1));
      c(j-i+1)=al*c(j-i+1)-be*conj(c(i));
      c(i)=top;
   end
   c(j)=-be/al;
end   
c=-c;
c(1)=1;
