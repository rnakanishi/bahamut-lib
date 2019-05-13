N = 32;
x = linspace(0, 1, N);
[X, Y, Z] = meshgrid(x, x, x);

figure
i = 5;
rotate3d on;
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
for i = 1:40
    v = dlmread(['../results/datatest/' int2str(i)]);
    v = reshape(v, N, N, N);
    v = permute(v, [2 3 1]);
    clf;
    % std(gcf,'Position', [10 10 900 600]);
    [face, vert] = isosurface(X, Y, Z, v, 0.0);
    % p = patch('faces', face, 'vertices', vert, 'EdgeColor', 'none');
    set(gcf, 'Position', [0 200 800 800])
    p = trisurf(face, vert(:, 1), vert(:, 2), vert(:, 3));
    light('Position', [2 1 5]);
    set(p, 'EdgeColor', 'none');
    set(p, 'FaceLighting', 'phong');
    title(int2str(i));
    view(75,30);
    axis([0 1 0 1 0 1], 'square');
    set(get(gca, 'XLabel'), 'String', 'Z axis');
    set(get(gca, 'YLabel'), 'String', 'X axis');
    set(get(gca, 'ZLabel'), 'String', 'Y axis');

    % axis equal tight
    pause;
    % saveas(gcf, ['animation3d/' int2str(i) '.png'])
end

% im = imread ('animation.pdf', 'Index', 'all');
% imwrite (im, 'animation.gif', 'DelayTime', .5)
