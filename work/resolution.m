function [R]=resolution(FH,F,Cd,Cm,rtmethod,p)
if nargin<5|isempty(rtmethod) rtmethod='PRT';end
np=length(Cm);

if nargin<6|isempty(p) p=0:np-1;end
fig=[221 222 223 224];
%fig=[223 224 223 224];
prow=np/2;
CdI=1./Cd;
CdI=diag(CdI);
CmI=1./Cm;
CmI=diag(CmI);
figure(10),
figure(11),

cd=1;
R=FH*cd*CdI*F;
figure(10);subplot(fig(1)),mesh(abs(R)),view(-15,30);
xlabel('columns'),ylabel('rows'),title('(a)');
figure(11),subplot(fig(1)),plot(p,abs(R(prow,:)));
xlabel('p'),ylabel('Amplitude'),title('(a)');


cm=1;cd=1;
R=inv(FH*cd*CdI*F+cm*CmI)*(FH*cd*CdI*F);
figure(10);subplot(fig(2)),mesh(abs(R)),view(-15,30);
xlabel('columns'),ylabel('rows'),title('(b)');
figure(11),subplot(fig(2)),plot(p,abs(R(prow,:)));
xlabel('p'),ylabel('Amplitude'),title('(b)');

cm=1e-5;cd=1;

R=inv(FH*cd*CdI*F+cm*CmI)*(FH*cd*CdI*F);
figure(10),subplot(fig(3)),mesh(abs(R)),view(-15,30);
xlabel('columns'),ylabel('rows'),title('(c)');
figure(11),subplot(fig(3)),plot(p,abs(R(prow,:)));
xlabel('p'),ylabel('Amplitude'),title('(c)');


if (rtmethod=='PRT')
cm=1;cd=1;
CmI(prow,prow)=1e-5;
%CmI(prow,prow)=1e-3;
R=inv(FH*cd*CdI*F+cm*CmI)*(FH*cd*CdI*F);
figure(10);subplot(fig(4)),mesh(abs(R)),view(-15,30);
xlabel('columns'),ylabel('rows'),title('(d)');
figure(11);subplot(fig(4)),plot(p,abs(R(prow,:)));
xlabel('p'),ylabel('Amplitude'),title('(d)');

end

