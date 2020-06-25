load xxc;
s=xxc(1:512,1);
ls = length(s);
 subplot(331);plot(s);title('s')

 [ca1,cd1] = dwt(s,'db1');
 % Perform one step reconstruction of    
	% ca1 and cd1.                          
                                        
   a1 = upcoef('a',ca1,'db1',1,ls);
   d1 = upcoef('d',cd1,'db1',1,ls);
   subplot(332);plot(ca1);title('ca1')
   subplot(333);plot(cd1);title('cd1')

subplot(334);plot(a1);title('a1')
subplot(335);plot(d1);title('d1')
subplot(336);plot(a1+d1);title('a1+d1')


a0 = idwt(ca1,cd1,'db1',ls);
[c,l] = wavedec(s,3,'db1');
% Extract approximation coefficients    
% at level 3, from wavelet decomposition
% structure [c,l].                      
                                        
ca3 = appcoef(c,l,'db1',3);  
subplot(337);plot(ca3);title('ca3')

% Extract detail coefficients at levels 
% 1, 2 and 3, from wavelet decomposition
% structure [c,l].                      
                                        
        cd3 = detcoef(c,l,3);           
        cd2 = detcoef(c,l,2);           
        cd1 = detcoef(c,l,1);          
subplot(338);plot(cd3);title('cd3')

% Reconstruct approximation at level 3,
% from the wavelet decomposition       
% structure [c,l].                     
                                       
a3 = wrcoef('a',c,l,'db1',3);  
subplot(339);plot(a3);title('a3')

% Reconstruct detail coefficients at 
% levels 1, 2 and 3, from the wavelet
% decomposition structure [c,l].     
                                     
        d3 = wrcoef('d',c,l,'db1',3);
        d2 = wrcoef('d',c,l,'db1',2);
        d1 = wrcoef('d',c,l,'db1',1);
% Reconstruct s from the wavelet
% decomposition structure [c,l].
                                
a0 = waverec(c,l,'db1');

figure,
subplot(331);plot(a0);title('a0')
subplot(332);plot(d1);title('d1')
subplot(333);plot(d2);title('d2')
subplot(334);plot(d3);title('d3')
subplot(334);plot(a0+d1+d2+d3);title('a0+d1+d2+d3')
