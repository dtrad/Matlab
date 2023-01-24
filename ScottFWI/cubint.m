function [ alpha ] = cubint( a , b , phia , phib , ga , gb )

    if a==b
        alpha=a;
    else
        d1 = ga + gb -3*(phia - phib)/(a - b);

        d2 = ( (d1^2) - ga*gb )^(1/2);

        alpha = b - (b - a)*(( gb + d2 - d1 )/( gb - ga + 2*d2));
    end

end

