clear;

folder = '~/git/bahamut-lib/results/particles/pls/';
timesteps = 100;
first = dlmread([folder 'ls' num2str(0)]);
N = sqrt(length(first(:, 1)));

x = linspace(-5, 5, N);
dx = x(2) - x(1);
[X, Y] = meshgrid(x, x);

az = 0;
el = 0;

fig = figure('position', [350, 150, 1500, 1500]);
% figure/
hold on;

% vis = [0 1 0 1];
vis = [-5 5 -5 5];

for t = 0:1:377
    % for t = [100 150 200]
    clf, hold on;
    % view(az, el);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    particles = dlmread([folder num2str(t)]);
    % escaped = dlmread([folder 'es' num2str(t)]);

    [row, col] = find(particles(:, 5) == -1);
    negatives = particles(row, :);
    [row, col] = find(particles(:, 5) == +1);
    positives = particles(row, :);

    hold on;

    % if (length(escaped) > 0)
    %  scatter(escaped(:, 1), escaped(:, 2), 50, 'g', '*');
    % end

    plot(negatives(1:3:end, 1), negatives(1:3:end, 2),'markersize', 0.5, 'b.');
    plot(positives(1:3:end, 1), positives(1:3:end, 2),'markersize', 0.5, 'r.');

    % quiver(particles(:, 1), particles(:, 2), particles(:, 2), -particles(:, 1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    levelset = dlmread([folder 'ls' num2str(t)]);
    values = reshape(levelset(:, 3), N, N)';
    levelx = reshape(levelset(:, 1), N, N)';
    hold on;
    levely = reshape(levelset(:, 2), N, N)';
    % contour(levelx, levely, values, -3 * dx:dx:0 * dx, 'c', 'linewidth', 2);
    % contour(levelx, levely, values, 0:dx:3 * dx, 'm', 'linewidth', 2);

    % subplot(121)
    hold on;
    contour(levelx, levely, values, [0, 0], 'k-', 'linewidth', 3);
    % subplot(122)
    % surf(levelx, levely, values);
    set(gca, 'xtick', levelx(1, :));
    set(gca, 'ytick', levely(:, 1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % grad = dlmread([folder 'g' num2str(t)]);
    % scale = 5;
    % gradx = reshape(grad(:, 3), N, N)';
    % grady = reshape(grad(:, 4), N, N)';
    % positionx = reshape(grad(:, 1), N, N)';
    % positiony = reshape(grad(:, 2), N, N)';
    % hold on;
    % quiver(positionx, positiony, scale * gradx, scale * grady, 'k');

    % [U, V] = gradient(values);
    % quiver(positionx, positiony, scale*U, scale*V, 'r');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % grad = dlmread([folder 'v' num2str(t)]);
    % scale = 2.5;
    % gradx = reshape(grad(:, 3), N, N)';
    % grady = reshape(grad(:, 4), N, N)';
    % positionx = reshape(grad(:, 1), N, N)';
    % positiony = reshape(grad(:, 2), N, N)';
    % hold on;
    % quiver(positionx, positiony, scale * gradx, scale * grady, 'y');

    % [U, V] = gradient(values);
    % quiver(positionx, positiony, scale*U, scale*V, 'r');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grid on;
    axis equal
    axis(vis);
    % axis([-5 5 -5 5]);
    title(t);

    axis off;
    pause(0);
    print(fig, ["../results/images/pls/" int2str(t) ".jpg"], "-S1150,1350")

    vis = axis();
    % pause;
    [az, el] = view();
end

levelset = dlmread([folder 'ls' num2str(0)]);
values = reshape(levelset(:, 3), N, N)';
levelx = reshape(levelset(:, 1), N, N)';
levely = reshape(levelset(:, 2), N, N)';
hold on;
contour(levelx, levely, values, [0, 0], 'k--', 'linewidth', 3);
