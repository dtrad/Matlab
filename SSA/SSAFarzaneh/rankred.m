%%d is the vector for each frequency, k is the number of linear events for rank reduction 
function [svdh]=rankred(d,k)
  h=hankel(d); 
    [U,S,V]=svd(h);
    S1=S;
    S1(k+1:end,k+1:end)=0;
    svdh=U*S1*V';
end
