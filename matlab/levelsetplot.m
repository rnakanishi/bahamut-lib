v = dlmread('levelset.txt');
N = sqrt(length(v));

x = linspace(-5, 5, N);
[X, Y] = meshgrid(x, x);

figure
i = 5;

% for i = 0:500
% z = dlmread(['../results/datatest/' int2str(i)]);
z = dlmread('levelset.txt');
z = reshape(z, N, N)'
clf;
contourf(X, Y, z, [-0 0]);
title(int2str(i));
axis([-5 5 -5 5], 'square');
set(gcf, 'Renderer', 'OpenGL');
set(gcf, 'Position', [0 0 800 601]);
% axis equal tight
% pause(1.0/60)
% print( ["animation/" int2str(i) ".png" ])
% saveas(gcf, ['animation/' int2str(i) '.png'])
% end

% im = imread ("animation.pdf", "Index", "all");
% imwrite (im, "animation.gif", "DelayTime", .5)
