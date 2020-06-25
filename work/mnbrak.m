function [a,b,c,fa,fb,fc]=mnbrak(a,b,fx,xpar,dpar,bpar,Apar,AApar,lambda)
test=0;
GOLD=1.618034;
GLIMIT=100.0;
TINY=1.0e-20;

if (test) 
  fa=feval(fx,a);
  fb=feval(fx,b);
else  
  fa=feval(fx,a,xpar,dpar,bpar,Apar,AApar,lambda,0);
  fb=feval(fx,b,xpar,dpar,bpar,Apar,AApar,lambda,0);
end


dum=0;

if (fb > fa) 
  [dum,a,b,dum]=shift4(dum,a,b,dum);
  [dum,fb,fa,dum]=shift4(dum,fb,fa,dum);
end

c=b+GOLD*(b-a);
if (test) fc=feval(fx,c);
else fc=feval(fx,c,xpar,dpar,bpar,Apar,AApar,lambda,0);end
while (fb > fc) 
  r=(b-a)*(fb-fc);
  q=(b-c)*(fb-fa);
  temp1=abs(max(abs(q-r),TINY))*sign(q-r);
  temp=sign(temp1);
  u=(b)-((b-c)*q-(b-a)*r)/(2.0*temp);
  ulim=(b)+GLIMIT*(c-b);
  if ((b-u)*(u-c) > 0.0) 
    if (test) fu=feval(fx,u);
    else fu=feval(fx,u,xpar,dpar,bpar,Apar,AApar,lambda,0);end
    if (fu < fc)
      a=b;
      b=u;
      fa=fb;
      fb=fu;
      return;
    elseif (fu > fb)
      c=u;
      fc=fu;
      return;
    end	
    u=(c)+GOLD*(c-b);
    if (test)     
      fu=feval(fx,u);
    else
      fu=feval(fx,u,xpar,dpar,bpar,Apar,AApar,lambda,0);
    end
  elseif ((c-u)*(u-ulim) > 0.0)
    if (test) 
      fu=feval(fx,u);
    else 
      fu=feval(fx,u,xpar,dpar,bpar,Apar,AApar,lambda,0);
    end
    if (fu < fc)
      [b,c,u]=shift4(b,c,u,c+GOLD*(c-b));
      if (test) futemp=feval(fx,u);
      else futemp=feval(fx,u,xpar,dpar,bpar,Apar,AApar,lambda,0);end
      [fb,fc,fu]=shift4(fb,fc,fu,futemp);
    end
  elseif ((u-ulim)*(ulim-c) >= 0.0)
    u=ulim;
    if (test) fu=feval(fx,u);
    else fu=feval(fx,u,xpar,dpar,bpar,Apar,AApar,lambda,0);end
  else
    u=c+GOLD*(c-b);
    if (test)       fu=feval(fx,u);
    else fu=feval(fx,u,xpar,dpar,bpar,Apar,AApar,lambda,0);end

  end	    
  [a,b,c,u]=shift4(a,b,c,u);
  [fa,fb,fc,fu]=shift4(fa,fb,fc,fu);
end


