function plotcost(JPF,J,JD,JP)
%subplot(221)
%sumiter=sum(JPF(1:60,:));semilogy(sumiter,'o');
%title('SUM of JPF from 1 to 60 Hz for each iteracion');
%xlabel('iteration');ylabel('sum JP'); 

subplot(211)
sumiter=sum(J(2:60,:));semilogy(sumiter,'o');
title('(a) Total Cost Function');
xlabel('iteration');ylabel('J'); 

subplot(223)
sumiter=sum(JD(2:60,:));semilogy(sumiter,'o');
title('(b) Data Misfit');
xlabel('iteration');ylabel('Jdata'); 

subplot(224)
sumiter=sum(JP(2:60,:));semilogy(sumiter,'o');
title('(c) Model Norm');
xlabel('iteration');ylabel('Jmodel'); 
