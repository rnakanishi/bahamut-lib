clear;
folder = '~/git/bahamut-lib/results/cip/3d/';
first = dlmread([folder num2str(0)]);
N = round(length(first)^(1/3));

timesteps = 73;

x = linspace(-5, 5, N);
[X, Y, Z] = meshgrid(x, x, x);

figure, hold on;
az = 0;
el = 90;

for count = 0:timesteps - 1
    clf, hold on;
    view(az, el);
    cip3d = dlmread([folder num2str(count)]);

    cip = reshape(cip3d, N, N, N);

    % anSurf = surf(X, Y, analytic);
    % weSurf = surf(X, Y, cip);
    hold on

    levels = 0;
    colors = colormap();
    colors = colors(1:floor(length(colors) / length(levels)):length(colors), :);

    for i = 1:length(levels)
        level = levels(i);
        isosurface(X, Y, Z, cip, level);
        % p = patch(isosurf);
        % isonormals(X, Y, Z, cip, p);
        % set(p, 'cdata', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', (1 / length(levels)) * i);
    end

    axis equal;
    axis([-5 5 -5 5 -5 5]);
    camlight

    % set(anSurf, 'cData', zeros(1, N * N));
    % set(weSurf, 'cData', ones(1, N * N));
    xlabel('Y');
    ylabel('X');
    zlabel('Z');
    title(num2str(count))

    pause();
    % pause;
    [az, el] = view();
end
