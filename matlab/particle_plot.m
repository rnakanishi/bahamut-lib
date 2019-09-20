clear;

folder = '~/git/bahamut-lib/results/particles/2d/';
timesteps = 100;
first = dlmread([folder 'ls' num2str(0)]);
N = sqrt(length(first(:, 1)));

x = linspace(-5, 5, N);
dx = x(2) - x(1);
[X, Y] = meshgrid(x, x);

az = 0;
el = 0;

figure('position', [350, 150, 1500, 1500]);
% figure/
hold on;

vis = [-5 5 -5 5];

for t = 0:100
    % for t = [0:34:1256 1256]
    clf, hold on;
    % view(az, el);
    particles = dlmread([folder num2str(t)]);
    escaped = dlmread([folder 'es' num2str(t)]);
    levelset = dlmread([folder 'ls' num2str(t)]);
    values = reshape(levelset(:, 3), N, N)';
    levelx = reshape(levelset(:, 1), N, N)';
    levely = reshape(levelset(:, 2), N, N)';

    [row, col] = find(particles(:, 5) == -1);
    negatives = particles(row, :);
    [row, col] = find(particles(:, 5) == +1);
    positives = particles(row, :);

    hold on;

    if (length(escaped) > 0)
        scatter(escaped(:, 1), escaped(:, 2), 50, 'g', '*');
    end

    scatter(positives(:, 1), positives(:, 2), 'b', '.');
    scatter(negatives(:, 1), negatives(:, 2), 'r', '.');

    % quiver(particles(:, 1), particles(:, 2), particles(:, 2), -particles(:, 1));

    contour(levelx, levely, values, -3 * dx:dx:0 * dx, 'm', 'linewidth', 2);
    contour(levelx, levely, values, 0:dx:3 * dx, 'c', 'linewidth', 2);
    contour(levelx, levely, values, [0, 0], 'k-', 'linewidth', 3);

    set(gca, 'xtick', levelx(1, :));
    set(gca, 'ytick', levely(:, 1));
    grid on;
    axis(vis);
    % axis([-5 5 -5 5]);
    title(t);
    pause;
    vis = axis();
    % pause;
    [az, el] = view();
end
