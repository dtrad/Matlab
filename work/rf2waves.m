function [R,a]=rf2waves(c,N)
% [R,a]=rf2waves(c,N)
c=-c;c(1)=1;R(1)=1;a(1)=1;V(1)=1;
for j=2:N
   a(j)=0;
   R(j)=c(j)*V(j-1);
   V(j)=V(j-1)*(1-c(j)*conj(c(j)));
   for i=2:j
      R(j)=R(j)-a(i)*R(j-i+1);
   end
   jh=(j+1)/2;
   for i=1:jh
      bot=a(j-i+1)-c(j)*conj(a(i));
      a(i)=a(i)-c(j)*conj(a(j-i+1));
      a(j-i+1)=bot;
   end   
end
R=-R;
R(1)=1;
