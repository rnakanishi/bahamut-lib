clear;
N = 40;
timesteps = 73;

x = linspace(-2, 2, N);
[X, Y, Z] = meshgrid(x, x, x);

figure, hold on;
az = 0;
el = 90;

for count = 0:timesteps - 1
    clf, hold on;
    view(az, el);
    weno3d = dlmread(['~/git/bahamut-lib/results/redistance/3d/' num2str(count)]);

    weno = reshape(weno3d, N, N, N);

    % anSurf = surf(X, Y, analytic);
    % weSurf = surf(X, Y, weno);
    hold on

    levels = -1:0.5:1;
    colors = colormap();
    colors = colors(1:floor(length(colors) / length(levels)):length(colors),:);

    for i = 1:length(levels)
      level = levels(i);
        isosurf = isosurface(X, Y, Z, weno, level);
        p = patch(isosurf);
        isonormals(X, Y, Z, weno, p);

        set(p, 'cdata', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', (1 / length(levels)) * i);
    end

    axis equal;
    axis([0 2 -2 2 -2 2]);
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
