function [V,DR]=radonf1(D,w,h,p)
% Daniel Trad-- 6-04-98
nh=length(h);
np=length(p);
DR=zeros(1,nh);
nf=length(w);
Qp=ones(1,np);
Qp(10)=0.001;Qp(19)=0.001;
ke=zeros(1,length(w));
index=0;
for f=2:nf-30,
   [FH,F]=radonmat2(w(f),h,p);
   v=V(f,:).';  %'
   dr=F*v;
   DR=[DR; dr.'];
end











