..it contiues from previous page. 



if (karmarkar)
  % Karmarkar interior point method
  method=['Log Barrier']
  iter_end=25;
  if (noisep==0)
    iter_end=15;
    PI_t=1e-5;
    DI_t=1e-5;
    DG_t=1e-5;
    gamma=1e-7;
    delta=1e-7;
    tol=1e-10;
    maxit=200;
  elseif (noisep==2.5)
    PI_t=1e-5;
    DI_t=1e-5;
    DG_t=1e-5;
    gamma=5e-2;
    delta=1e-1;
    tol=1e-6;
    maxit=70;
  elseif (noisep==5)
    iter_end=15;
    PI_t=1e-6;
    DI_t=1e-6;
    DG_t=1e-6;
    gamma=0.001e-1;
    delta=6e-1;
    tol=1e-6;
    maxit=150;
  end  
    
  figure(4);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1),title(['Log Barrier ']);
  hold on
  tic
  x_sol=karmarkar_interface(d,AE,iter_end,PI_t,DI_t,DG_t,gamma, ...
			     delta,tol,maxit);
end

if (nlcg1)
  method=['Nonlinear CG']
  fx='l2l1cost';
  lambda=1e-1;
  tol=1e-5;
  itercg=1000;%*length([x;x2]);
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  tic
  [x_sol,delnew,niter]=nlcgproject4(fx,sparse(AE),d,tol,x0,Wm,lambda,itercg,0,xt);
  
end

if (nlcg2)
  method=['Nonlinear CG']
  fx='l2l1cost'
  lambda=1e-2;
  itercg=3*length([x;x2]);
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  [x_sol,fy,normdeltay,normgrad]=nlcgproject5(fx,x0,d,sparse(AE),lambda,itercg,tol); % NL CG
  
end
if (nlcg3)
  method=['Nonlinear CG']
  fx='l2l1costparam';
  if (noisep==0) lambda=1e-0;tol=1e-5;itercg=100
  elseif (noisep==2.5) lambda=5e-0;tol=1e-2;itercg=200
  elseif (noisep==5) lambda=5e-2;tol=1e-5;itercg=300
  end
  %lambda=1e-0;
  
  Wm=eye(nmodel);
  figure(5);
  subplot(1,1,1)
  plot(1:length([x;x2]),[x;x2],'x'),axis(axis1), 
  title(['Non-linear CG']);
  %hold on
  tic
  [x_sol,delnew,niter]=nlcgproject6(fx,sparse(AE),d,tol,x0,Wm,lambda,itercg,1,xt);
    
end

  
text1=sprintf('Noise=%4.1f %% , required time %5.2f s\n',noisep, toc);
dp=AE*x_sol;
bp=A*x_sol(1:length(x));
ep=AH*x_sol(length(x)+1:length(xt));

hold off  
plot(1:length(x_sol),x_sol,'r'),axis(axis1), 

figure(10)
subplot(311)
plot(1:nmodel,x_sol,1:nmodel,[x;x2],'.'),     
axis(axis1),title([method,' : ',text1])
subplot(312)
plot(1:length(b),b,'.',1:length(b),bp),axis(axis2);  
title('signal and predicted signal')    
subplot(313)
plot(1:length(e),e,'.',1:length(e),ep),axis(axis2);  
title('noise and predicted noise')    


return;
