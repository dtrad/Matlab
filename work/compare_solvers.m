N=30:100;
C=N.^3/6;
L=N.^2;
L2=N.*log(N).^2;
figure,
plot(N,C,'o',N,L,'+');
title('Cholesky (o), Levinson (+)');
ylabel('Number of operations');
xlabel('Size of the system');

