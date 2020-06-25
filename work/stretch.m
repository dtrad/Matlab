        function [d]=t_2(dt,dts,d,index)
        
% 		 function [d]=t_2(dt,dts,d,index)
%       Index.eq. 1 then stretch the time axis     t-->t'=t^2
%       Index.eq.-1 then unstretch the time axis
 
 
%       when INDEX=1 this is the IN/OUT
 
%       Input parameters:
%
%         d(nt,nh)    - data before tranformation
%               ds    - original sampling rate (sec)
%
%        Output parameters
%
%         d(nts,nh)  - data after tranformation
%         dts        - sampling rate in the stretched axis (sec**2)
%
%        Note:
%
%        When index=-1 the IN/OUT are interchanged

%        M.D.Sacchi, Physics, UofA. 
 

[nt nh]=size(d);

ih=1:nh;
    for i=1:nt;
        if(index.eq. 1) t=(sqrt(i*dts))/dt;end
        if(index.eq.-1) t=((i*dt)^2)/dts;end
        if(t>1) 
        		i1=fix(t);
        		i2=i1+1;
        aux(i)=d(i1,ih)+(d(i2,ih)-d(i1,ih))*(t-i1)/(i2-i1);
        else
        aux(i)=0.d0;
        end
     end
i=1:nt;d(i,ih)=aux(i);   
    
        
% --------------------------------------------------------------