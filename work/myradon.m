function [R,xp]=myradon(I)
figure,
subplot(221);image(I);colorbar;title('(a)');

        theta = 0:3:180;
        [R,xp] = radon(I,theta);
        subplot(222);
        imshow(theta,xp,R*10,[],'n')
        title('(b)');
        xlabel('\theta (degrees)')
        ylabel('x''')
        %colormap(hot),
        colorbar
        figure(gcf)