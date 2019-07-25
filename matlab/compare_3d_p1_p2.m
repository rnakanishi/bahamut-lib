v = dlmread('../results/datatest/0');
N = round(length(v)^(1/3));
v = dlmread('../results/data3d-40/datatest/0');
N2 = round(length(v)^(1/3));

x = linspace(0, 1, N);
[X, Y, Z] = meshgrid(x, x, x);
x = linspace(0, 1, N2);
[X2, Y2, Z2] = meshgrid(x, x, x);

figure
i = 5;
rotate3d on;
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
for i = 0:132
    [az, el] = view();
    v = dlmread(['../results/datatest/' int2str(i)]);
    v = reshape(v, N, N, N);
    v = permute(v, [2 3 1]);
    clf;
    % std(gcf,'Position', [10 10 900 600]);
    [face, vert] = isosurface(X, Y, Z, v, 0.0);
    % p = patch('faces', face, 'vertices', ver55t, 'EdgeColor', 'none');
    subplot(121)
    set(gcf, 'Position', [0 200 900 900])
    p = trisurf(face, vert(:, 1), vert(:, 2), vert(:, 3));
    light('Position', [2 1 5]);
    set(p, 'EdgeColor', 'none');
    set(p, 'FaceLighting', 'phong');
    title(int2str(i));
    view(az, el);
    axis([0 1 0 1 0 1], 'square');
    set(get(gca, 'XLabel'), 'String', 'Z axis');
    set(get(gca, 'YLabel'), 'String', 'X axis');
    set(get(gca, 'ZLabel'), 'String', 'Y axis');


    v = dlmread(['../results/data3d-40/datatest/' int2str(i)]);
    v = reshape(v, N2, N2, N2);
    v = permute(v, [2 3 1]);
    subplot(122);
    % std(gcf,'Position', [10 10 900 600]);
    [face, vert] = isosurface(X2, Y2, Z2, v, 0.0);
    % p = patch('faces', face, 'vertices', ver55t, 'EdgeColor', 'none');
    set(gcf, 'Position', [0 200 900 900])
    p = trisurf(face, vert(:, 1), vert(:, 2), vert(:, 3));
    light('Position', [2 1 5]);
    set(p, 'EdgeColor', 'none');
    set(p, 'FaceLighting', 'phong');
    title(int2str(i));
    view(az, el);
    axis([0 1 0 1 0 1], 'square');
    set(get(gca, 'XLabel'), 'String', 'Z axis');
    set(get(gca, 'YLabel'), 'String', 'X axis');
    set(get(gca, 'ZLabel'), 'String', 'Y axis');
    % axis equal tight
    pause(0.05);
    % saveas(gcf, ['animation3d/' int2str(i) '.png'])
end

% im = imread ('animation.pdf', 'Index', 'all');
% imwrite (im, 'animation.gif', 'DelayTime', .5)