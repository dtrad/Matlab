function [x]=karmarkar_interface(b,A,iter_end,PI_t,DI_t,DG_t,gamma,delta)
AA=[A -A];
[mm]=karmarkar(b,AA,iter_end,PI_t,DI_t,DG_t,gamma,delta);
n=length(mm);
x=mm(1:n/2)-mm(n/2+1:end);