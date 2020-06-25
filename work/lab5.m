% Lab # 5

nt=100;
w1=pi/4;
w2=pi/2;
dt=1;
t=0:nt-1;
t=t*dt;
x=sin(w1*t)+sin(w2*t);
rhob=1;
rhoa=1.2;

f=-nt/2:nt/2-1;
f=f*2/nt*1/(2*dt);
f=f*2*pi;

figure(1)
subplot(311),plot(t,x);title('x')
subplot(312),plot(f,fftshift(abs(fft(x))));
subplot(313),plot(f,fftshift(angle(fft(x))));

rhob=1
answer='y'
while answer=='y' 
    B=[1 -2/rhob*cos(pi/4) 1/rhob^2];
    Bz=[B zeros(1,length(x)-length(B))];


    figure(2)
    subplot(211),plot(f,fftshift(abs(fft(Bz))));title('FIR filter')
    subplot(212),plot(f,fftshift(angle(fft(Bz))));
    figure(gcf)
    answer=input('another try?','s');
    if (answer=='y')
       rhob=input('rhob=?'),
    end
end


%y=conv(x,B);
y=filter(B,1,x);
y=y(1:nt);

figure(3)
subplot(311),plot(t,y);title('y=conv(b,x)')
subplot(312),plot(f,fftshift(abs(fft(y))));
subplot(313),plot(f,fftshift(angle(fft(y))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IIR
rhoa=1
answer='y'
while answer=='y' 

     A=[1 -2/rhoa*cos(pi/4) 1/rhoa^2];
     Az=[A zeros(1,length(x)-length(A))];

     figure(4)
     subplot(211),plot(f,fftshift(abs(fft(Bz)./fft(Az))));title('B/A')
     subplot(212),plot(f,fftshift(angle(fft(Bz)./fft(Az))));title('phase')

     y=filter(B,A,x);
     y=y(1:nt);
     figure(gcf)

     figure(5)
     subplot(311),plot(t,y);title('y=conv(b,x)')
     subplot(312),plot(f,fftshift(abs(fft(y))));
     subplot(313),plot(f,fftshift(angle(fft(y))));



     figure(gcf)
     answer=input('another try?','s');
     if (answer=='y')
       rhoa=input('rhoa=?'),
     end
end













