function [time] = dsqr(tau,off1,off2,vel)
[time]=sqrt(tau*tau/4+off1*off1/vel/vel)+sqrt(tau*tau/4+off2*off2/vel/vel);
return; 