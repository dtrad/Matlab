clear
close all;
figure
Z = peaks;
surf(Z)
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';

vidObj= VideoWriter('peaks.avi');
open(vidObj);

loops = 40;
F(loops) = struct('cdata',[],'colormap',[]);
for j = 1:loops
    X = sin(j*pi/10)*Z;
    surf(X,Z)
    drawnow
    currFrame=getframe;
    writeVideo(vidObj,currFrame);
    F(j) = getframe;
end
close(vidObj);
