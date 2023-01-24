%%svdh is the rank reduced hankel matrix, mid is the dimendions of the hankel matrix, lent is the length of the vector
function [s]=antidiag(svdh,mid,lent)
  for i=1:mid
    if i+mid-1<=lent
      s(i)=svdh(1,i);
    if i>=2
      s(i+mid-1)=svdh(i,mid);
    end
    end
  end 
end