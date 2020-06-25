function [f,F,y2,y3]=shapingfilter2d(y,x,epsilon)
% Use Wienner to calculate the shaping filter
% that transform d into dout
% d data (input to the filter)
% dout desired output
% epsilon regularization

% y: freq domain filtering
% y2: time domain filtering


%load sudata.mat

if (nargin<3) epsilon=1;end
if (size(x)~=size(y)) 
  echo 'size of d and dout has to be the same'
  return
end  

[n1 n2]=size(x);

X=fft2(x);
Y=fft2(y);

F=(conj(Y).*X./(max(epsilon,(conj(Y).*Y))));
Y2=Y.*F;

Y2=duplic2d(Y2);
y2=real(ifft2(Y2));

F=duplic2d(F);
f=real(fftshift(ifft2(F)));


% Time domain 2d deconvolution
y3=conv2(y,f);

n12=n1/2;
n22=n2/2;

y3=y3(n12:n12+n1-1,n22:n22+n2-1);


return;



