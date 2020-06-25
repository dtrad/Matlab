% In this script I test a 3-dimensional Curvelet

% Parameters for the curvelets
dlevs = [1,1,2,2,3];
pfilt = 'db8';
dfilt = 'pkva';

% Make a curvelet vector
N        =  256;
ndim     =  2;
cx       =  zeros(N);
cx       =  pdfbdec(cx,pfilt,dfilt,dlevs);

[vcy, s] =    pdfb2vec(cx);
fssize   =    prod(s(end, :));	
vcy(fssize+1000)  = 1;

cy       =    vec2pdfb(vcy, s);

y        =    pdfbrec(cy,pfilt, dfilt);


imagesc(y)
%y=zeros(size(y));
%y(128,128)=1;
%w=ricker(35,0.004);


%yy=conv(y(:,128),w);
%yy=yy(1:N);
%y(:,128)=yy(:);




writesudata('/home/dtrad/work/onecurve.bin',y)

eval(['!cd /home/dtrad/work;Ximageonecurvelet;cd ../matlab']);

[yadj]=readsudata('/home/dtrad/work/onecurve.adj.bin',N,N);




