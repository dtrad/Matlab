 function plotll(LLO,f) 
   if (f==3)
   figure
   subplot(311),plot(svd(LLO));title('(a)'), 
   elseif (f==50)
   subplot(312),plot(svd(LLO));title('(b)'), 
   elseif (f==100)
   subplot(313),plot(svd(LLO));title('(c)'),      
   end

