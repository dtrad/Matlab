function mycomplexplot(c,xaxis)
if nargin < 2 xaxis=1:length(c);end
figure,
subplot(411);plot(xaxis,real(c));title('real');ylabel('amplitude');xlabel('p(s/m)')
subplot(412);plot(xaxis,imag(c));title('imag');ylabel('amplitude');xlabel('p(s/m)')

%subplot(413);plot(xaxis,exp(-(imag(c))));title('exp(-(imag(c)))')
%subplot(414);plot(xaxis,exp(-abs(imag(c))));title('exp(-abs(imag(c)))')
