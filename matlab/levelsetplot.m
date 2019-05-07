x = linspace(0, 1, 256);
[X, Y] = meshgrid(x, x);

figure
i = 5;

for i = 0:300
    z = dlmread(['data/' int2str(i)]);
    clf;
    contourf(X, Y, z, [-0 0]);
    title(int2str(i));
    axis([0 1 0 1],'square');
    set(gcf, 'Renderer', 'OpenGL');
    set(gcf,'Position',[0 0 800 601]);
    % axis equal tight
    pause(0.15)
    % print( ["animation/" int2str(i) ".png" ])
    saveas(gcf, ['animation/' int2str(i) '.png'])
end

% im = imread ("animation.pdf", "Index", "all");
% imwrite (im, "animation.gif", "DelayTime", .5)
