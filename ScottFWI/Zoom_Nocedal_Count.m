function [alphastar,eval] = Zoom_Nocedal_Count( x0 , descentd , c1 , c2 , fh , phi0 , g0 , alphalo , alphahi , eval )
%%
i=1;

found=0;

[philo,glo] = fh(x0 + alphalo*descentd);
eval = eval + 1 ;

glo = glo'*descentd;

[phihi,ghi] = fh(x0 + alphahi*descentd);
eval = eval + 1 ;

ghi = ghi'*descentd;

save initial alphalo  alphahi  philo  phihi  glo  ghi descentd x0
    
    while i <= 10
        
        
        if alphalo < alphahi
            
            alpha=cubint( alphalo , alphahi , philo , phihi , glo , ghi );
            
            if imag(alpha) ~= 0
                save cubinterr alphalo alphahi philo phihi glo ghi 
                break
            end
            
        elseif alphalo == alphahi
            
            alpha=alphalo; 
        
        else
            
            alpha=cubint( alphahi , alphalo , phihi , philo , ghi , glo );
            
            if imag(alpha) ~= 0
                save cubinterr alphalo alphahi philo phihi glo ghi 
                break
            end
            
        end
        
        
        
        if isnan(alpha)
            
            save bad alphalo  alphahi  philo  phihi  glo  ghi descentd x0 i
            
            break
            
        end
        
        [phi,g] = fh(x0 + alpha*descentd);
        eval = eval + 1 ;
    
        g = g'*descentd;

        if phi > phi0 + c1*alpha*g0 || phi >= philo

            alphahi = alpha;
            
            phihi = phi;
            
            ghi = g;

        else

            if abs(g) <= -c2*g0

                alphastar = alpha;

                found=1;

                break

            end

            if g*(alphahi - alphalo) >= 0

                alphahi = alphalo;

                phihi = philo;

                ghi = glo;

            end

            alphalo = alpha;

            philo = phi;

            glo = g;
        
        end
        
        i=i+1;
        
        if alphalo == alphahi
            
            alphastar=alphalo;
            
            break
        end
    
    end
    
    if ~found
    
        if alphalo < alphahi
            
            alphastar = cubint( alphalo , alphahi , philo , phihi , glo , ghi );
            
            if isnan(alphastar)
                
                save bad alphalo  alphahi  philo  phihi  glo  ghi descentd x0 i
      
            end
            
        else
            
            alphastar = cubint( alphahi , alphalo , phihi , philo , ghi , glo );
            
            if isnan(alphastar)
                
                save bad alphalo  alphahi  philo  phihi  glo  ghi descentd x0 i
                
            end
    
        end
    end

end

