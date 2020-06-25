% Image coding.                       
nbcol = size(xx,1);          
family='db4';

% Perform one step decomposition.                
[ca1,ch1,cv1,cd1] = dwt2(xx,family);       
figure,wigb([ca1,ch1;cv1,cd1]);

% Perform second step decomposition :       
% decompose approx. cfs of level 1.         
[ca2,ch2,cv2,cd2] = dwt2(ca1,family);
                                            
figure,wigb([ca2,ch2;cv2,cd2]);

% Invert directly decomposition of X              
% using coefficients at level 1.                  
                                                  
a0 = idwt2(ca1,ch1,cv1,cd1,family,size(xx));
figure,wigb(a0);

% Perform decomposition at level 2  
% of X using db1.                   
                                    
[c,s] = wavedec2(xx,2,family);
% Extract coefficients at level 2,           
% from wavelet decomposition structure [c,s].
                                             
        ca2 = appcoef2(c,s,family,2);         
        ch2 = detcoef2('h',c,s,2);           
        cv2 = detcoef2('v',c,s,2);           
        cd2 = detcoef2('d',c,s,2);           
        % Extract coefficients at level 1,           
% from wavelet decomposition structure [c,s].
                                             
        ca1 = appcoef2(c,s,family,1);         
        ch1 = detcoef2('h',c,s,1);           
        cv1 = detcoef2('v',c,s,1);           
        cd1 = detcoef2('d',c,s,1);        
        
        
% Reconstruct approximation at level 2,
% from the wavelet decomposition       
% structure [c,s].                     

a2 = wrcoef2('a',c,s,family,2); 
% Reconstruct detail at level 2,      
% from the wavelet decomposition      
% structure [c,s].                    
        h2 = wrcoef2('h',c,s,family,2);
        v2 = wrcoef2('v',c,s,family,2);
        d2 = wrcoef2('d',c,s,family,2);
        % One step reconstruction of wavelet
% decomposition structure [c,s].    
[c,s] = upwlev2(c,s,family); 
% Reconstruct approximation and details
% at level 1, from coefficients.       
%                                      
% step 1 : extract coefficients        
% decomposition structure [c,s].       
%                                      
% step 2 : reconstruct.     
 siz = s(size(s,1),:);              
    ca1 = appcoef2(c,s,family,1);       
    ch1 = detcoef2('h',c,s,1);         
    cv1 = detcoef2('v',c,s,1);         
    cd1 = detcoef2('d',c,s,1);         
    a1  = upcoef2('a',ca1,family,1,siz);
    hd1 = upcoef2('h',ch1,family,1,siz);
    vd1 = upcoef2('v',cv1,family,1,siz);
    dd1 = upcoef2('d',cd1,family,1,siz);
    % Reconstruct X from the wavelet 
% decomposition structure [c,s]. 
                                 
        a0 = waverec2(c,s,family);