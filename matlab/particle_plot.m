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

% vis = [0 1 0 1];
vis = [-5 5 -5 5];

for t = 70:1:630
    % for t = [100 150 200]
    clf, hold on;
    % view(az, el);

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

    scatter(positives(:, 1), positives(:, 2), 1.5, 'b', '.');
    scatter(negatives(:, 1), negatives(:, 2), 1.5, 'r', '.');

    % quiver(particles(:, 1), particles(:, 2), particles(:, 2), -particles(:, 1));

    levelset = dlmread([folder 'ls' num2str(t)]);
    values = reshape(levelset(:, 3), N, N)';
    levelx = reshape(levelset(:, 1), N, N)';
    levely = reshape(levelset(:, 2), N, N)';
    contour(levelx, levely, values, -3 * dx:dx:0 * dx, 'm', 'linewidth', 2);
    contour(levelx, levely, values, 0:dx:3 * dx, 'c', 'linewidth', 2);

    % subplot(121)
    contour(levelx, levely, values, [0, 0], 'k-', 'linewidth', 3);
    % subplot(122)
    % surf(levelx, levely, values);

    set(gca, 'xtick', levelx(1, :));
    set(gca, 'ytick', levely(:, 1));
    grid on;
    axis equal
    axis(vis);
    % axis([-5 5 -5 5]);
    title(t);

    pause();
    % print(["../results/particleLevelset/weno" int2str(t) ".png"])

    vis = axis();
    % pause;
    [az, el] = view();
end
