load ozdataresult;
figure(1);clf
subplot(211);dd=d(1:512,:);dd(:,24:25)=dd(:,24:25)*0.1;wigb(dd,10,h,t(1:512),1);
title('Offset(m)');ylabel('Time(s)');
subplot(212);dd=d2(1:512,:);dd(:,24:25)=dd(:,24:25)*0.1;wigb(dd,10,h,t(1:512),1);
title('Offset(m)');ylabel('Time(s)');