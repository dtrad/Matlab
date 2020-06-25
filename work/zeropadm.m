function [outm]=zeropadm(in,lx)
% [out]=padzeros(in,lx)
% pad zeros on IN till length lx
% Daniel Trad- UBC- 23/07/98
in=seis_shape(in);
[NT,NX]=size(in);
for ii=1:NX
   clear out;
   out=in(:,ii);
   out(:);
   ll=length(out);
   out=[out;zeros(lx-ll,1)];
   if ii~=1 
      outm=[outm out];
     else
      outm=out;
   end   
end



