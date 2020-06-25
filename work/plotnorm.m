function plotnorm(sigmap,J,JP,JD)
% Plot misfits and norms
% Daniel Trad- 6-04-98
figure,
subplot(221),loglog(sigmap,real(J),'o');title('sigmap-J');
xlabel('sigmap');ylabel('Total norm');
subplot(222),loglog(sigmap,real(JD),'o');title('sigmap-JD');
xlabel('sigmap');ylabel('Data misfit');
subplot(223),loglog(sigmap,real(JP),'o');title('sigmap-JP');
xlabel('sigmap');ylabel('Model norm');
subplot(224),loglog(real(JP),real(JD),'o');title('JP-JD');
xlabel('model norm');ylabel('Data misfit');
