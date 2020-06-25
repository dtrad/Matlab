function c=vector2matrix(c,dt1,dt,NF)
% c=vector2matrix(c,dt1,dt,NF)
   NP=length(c)
   ctemp=zeros(NF,NP);
   for pp=1:NP,
      dt1i=t2index(dt1(pp),dt);
      ctemp(dt1i,pp)=c(pp);
   end
   c=ctemp;
