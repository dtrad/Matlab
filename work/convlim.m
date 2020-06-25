function ab=convlim(a,b,c,ab) 
% Convolution between a and b modifies the c first elements
% of ab;
% ab=convlim(a,b,c,ab)
% Daniel Trad- UBC
mm=length(a)+length(b)-1;
if nargin < 3 c=mm;end
if nargin < 4 ab=zeros(c,1);end;

temp=conv(a,b);
if c>mm temp(mm+1:c)=zeros(c-mm,1);end;
ab(1:c)=temp(1:c);
