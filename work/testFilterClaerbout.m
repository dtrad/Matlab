%program to test triangular filter
close all
if (1==1)
            
load xtest
t=1:length(x);
figure(1);
plot(t,boxcarClaerbout(5,length(x),x),t,x)

figure(2);
plot(t,triangleClaerbout(5,length(x),x),t,x*1)

dt=0.004;
w=freqaxis(dt,length(x));
figure(3);
plot(w,fftshift(abs(fft(boxcarClaerbout(5,length(x),x)))),w,fftshift(abs(fft(x))))

figure(4);
plot(w,fftshift(abs(fft(triangleClaerbout(5,length(x),x)))),w,fftshift(abs(fft(x))))

end

%test with myfilter
figure(5);t=1:length(x);plot(t,triangleFilter(9,length(x),x),t,x*1)
figure(6);t=1:length(x);plot(w,fftshift(abs(fft(triangleFilter(5,length(x),x)))),w,fftshift(abs(fft(x))))