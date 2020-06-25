load stransfex
[S,s]=stransf0(x);
figure(1);subplot(311);plot(x);
figure(1);subplot(312);imagesc(abs(S));
figure(1);subplot(313);plot(real(s));