 function plotll(LLO,f) 
   if (f==3)
   figure
   subplot(321),mesh(real(LLO));title('(a)'),
   subplot(322),mesh(imag(LLO));title('(b)'), 
   elseif (f==50)
   subplot(323),mesh(real(LLO));title('(c)'),
   subplot(324),mesh(imag(LLO));title('(d)'), 
   elseif (f==100)
   subplot(325),mesh(real(LLO));title('(e)'),
   subplot(326),mesh(imag(LLO));title('(f)'),      
   end

