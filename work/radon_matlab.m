 I = zeros(100,100);
 I(25:75, 25:75) = 1;
 theta = 0:180;
 [R,xp] = radon(I,theta);
 imshow(theta,xp,R,[],'n')
 xlabel('\theta (degrees)')
 ylabel('x''')
 colormap(hot), colorbar
