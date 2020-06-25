function [img,tcoef]=recons(y,coef,btree,D,perc)

if (1)
    % Thresholding 
    % Define the threshold such that only 100 coefficients remain
    [N1 M1]=size(coef);
    total=N1*M1;
    kept=round((total)*(perc/100));
    if (kept>total) 
        kept=total;
    end
    message=sprintf('number of coefficients to keep=%d from %d or %f percent',...
        kept,total,kept/total*100)



    ss=sort((abs(coef(:))));q=ss(length(ss)-kept+1);
    tcoef=coef;tcoef(find(abs(coef)<q))=0;
    figure(3);imagesc(tcoef);title('thresholded ridgelet coefficients');colorbar;
    figure(gcf)
end

img = IPT2_RPkt(btree,tcoef,D);


figure;imagesc(y);title('y');
figure;imagesc(img);title('img')
figure;imagesc(y-img);title('y-img');

return;