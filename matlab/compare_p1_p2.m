v = dlmread('../results/datatest/0');
N = (length(v));

x = linspace(0, 1, N);
[X, Y] = meshgrid(x, x);

figure
i = 5;

for i = 0:250
    clf;
    z = dlmread(['../results/datatest/' int2str(i)]);
    subplot(121)
    contourf(X, Y, z, [-0 0]);
    title(['P1 ' int2str(i)]);
    axis([0 1 0 1], 'square');
    set(gcf, 'Renderer', 'OpenGL');
    set(gcf, 'Position', [0 0 800 601]);
    % axis equal tight

    z = dlmread(['../results/2d-p2/' int2str(i)]);
    subplot(122)
    contourf(X, Y, z, [-0 0]);
    title(['P2 ' int2str(i)]);
    axis([0 1 0 1], 'square');
    set(gcf, 'Renderer', 'OpenGL');
    set(gcf, 'Position', [0 0 800 601]);
    pause(1.0/60)
    % print( ["animation/" int2str(i) ".png" ])
    % saveas(gcf, ['animation/' int2str(i) '.png'])
end

% im = imread ("animation.pdf", "Index", "all");
% imwrite (im, "animation.gif", "DelayTime", .5)
