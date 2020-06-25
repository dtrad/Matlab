t=-128:127;t=t/100;
figure(1)
subplot(211),plot(t,t.^2,'o',t,eps+abs(t),'+',t,log(10+t.^2),t,log(2+t.^2))
xlabel('model');ylabel('Model Norm');title('(a)');
subplot(212),plot(t,2*t,'o',t,eps+sign(t),'+',t,2./(10+t.^2).*t,t,2./(2+t.^2).*t)
xlabel('model');ylabel('Gradient');title('(b)');
