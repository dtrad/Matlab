function [dataout]=correlmatrix(datain1, datain2)

[n1 n2]=size(datain1);
[n12 n22]=size(datain2);

if (n1~=n12 | n2 ~= n22 ) return; end


n1h=round(n1/2);
n2h=round(n2/2);

datain1z=[ zeros(n1h,2*n2); zeros(n1,n2h) datain1 ...
	  zeros(n1,n2h); zeros(n1h,2*n2)];

datain2z=[ zeros(n1h,2*n2); zeros(n1,n2h) datain2 ...
	  zeros(n1,n2h); zeros(n1h,2*n2)];

close all;


DATAKF1=fft2(datain1z);
DATAKF2=fft2(datain2z);

dataout=ifft2(DATAKF1.*conj(DATAKF2));

% Shifting

A=dataout;
dataout=[A(n1:end,:);A(1:n1-1,:)];

A=dataout;
dataout=[A(:,n2:end),A(:,1:n2-1)];

A=dataout;
dataout=A(n1h:n1h+n1,:); 

A=dataout;
dataout=A(:,n2h:n2h+n2); 

