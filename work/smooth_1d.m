      function [out]=smooth_1d(q, nl)
%       simple smoother
%       it does not change the borders.
%       q 1d vector size ns
%       qtemp 1d vector size ns to store old samples
%       ns size of q
%       nl length of smoother (odd number)

        ns=length(q);      
        if (ns <= nl)
            return
        end
            
        nl2 = (nl-1)/2;
        sum = 0;
        for i = 1:nl
           sum = sum + q(i);
        end

%       store original data
        qtemp(1:ns) = q(1:ns);
% print*,'smoothing ns = ',ns,' samples with nl = ',nl
        for i = nl2+2: ns - nl2
           sum  = sum - qtemp(i-nl2-1) + qtemp(i+nl2);
           q(i) = sum / nl;
        end
        out = q;
        return
         