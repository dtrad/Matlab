
function [alphastar,phi,eval,path] = Linesearch_Nocedal_Full( x0 , descentd , c1 , c2 , fh , alpha0 )
% NowScott: Scott, it's not working, I'm getting complex step lengths, are you
% sure you didn't screw up??
% PastScott: Is your descent direction actually an ascent direction?
% NowScott: ...
% This blurb has been added so that I don't have to check this code another
% 80 times
alpha=alpha0;

alphaprev=0;

[phi0,g0] = fh(x0);

phiprev = phi0;
% save check
g0 = g0'*descentd;

if g0 < 0
    
%     phis = fh(x0+1e-9*alpha0*descentd);
%     
%     silly = (phis - phi0)/(1e-9*alpha0);
%     
%     fprintf('1D derivative is supposed to be %5.5e \n', g0);
%     fprintf('But is actually %5.5e \n', silly);
    
    i=1;
    
    found = 0;
    
    alphastar = 0;
    
    eval=0;
    
    % save temp fh
    
    while i<=10
        
        [phi,g] = fh(x0 + alpha*descentd);
        eval = eval + 1 ;
        
        whileit=0;
        
        while phi == Inf
            whileit=whileit+1;
            alpha=alpha/2;
            [phi,g] = fh(x0 + alpha*descentd);
            eval = eval + 1 ;
            if whileit >= 10
                phi=0;
            end
        end
        
        if whileit >= 10
            break
        end
        
        
        g = g'*descentd;
        
        if phi > phi0 + c1*alpha*g0 || (phi >= phiprev && i > 1)
            
            [alphastar,eval] = Zoom_Nocedal_Count( x0 , descentd , c1 , c2 , fh , phi0 , g0 , alphaprev , alpha , eval );
            %         save temp x0 descentd c1 c2 fh phi0 g0 alphaprev alpha
            found = 1;
            
            path=1;
            
            break
            
        end
        
        if abs(g) <= -c2*g0
            
            alphastar = alpha;
            
            found = 1;
            
            path=2;
            
            break
            
        end
        
        if g >= 0
            
            [alphastar,eval] = Zoom_Nocedal_Count( x0 , descentd , c1 , c2 , fh , phi0 , g0 , alpha , alphaprev , eval );
            
            found = 1;
            
            path=3;
            
            break
            
        end
        
        alphaprev=alpha;
        
        alpha=2*alpha;
        
        phiprev=phi;
        
        i=i+1;
        
        if imag(alphastar)~=0
            save imag_err path alpha
            break
        end
        
    end
    
    if ~found
        alphastar=0;
        path=4;
    end
else
    alphastar = 0;
    eval = 1;
    phi = phi0;
    path = 5;
end

neg_desc = descentd < 0;
ratio = x0./descentd;
alpha_zero = min(abs(ratio(find(neg_desc~=0))));
alpha_max = 0.9*alpha_zero;
if alphastar > alpha_max
    alphastar = alpha_max;
end

[phi] = fh(x0 + alphastar*descentd);