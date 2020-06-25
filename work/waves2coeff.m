function [c]=waves2coeff(a,n)
c(1)=-1;R(1)=1;a(1)=1;V(1)=1;
i=1:n;c(i)=a(i);
for k=1:n
   j=n-k+2;   
   al=1/(1-c(j)*conj(c(j)));
   be=c(j)*al;
   jh=fix((j+1)/2);
   for i=1:jh
      top=al*c(i)-be*conj(c(j-i+1));
      c(j-i+1)=al*c(j-i+1)-be*conj(c(i));
      c(i)=top;
   end
   c(j)=-be/al;
end   
  