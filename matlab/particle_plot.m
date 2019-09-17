clear;

folder = '~/git/bahamut-lib/results/particles/2d/';
timesteps = 100;
first = dlmread([folder 'ls' num2str(0)]);
N = sqrt(length(first));

x = linspace(-5, 5, N);
[X, Y] = meshgrid(x, x);
dx = x(2) - x(1);

az = 0;
el = 0;

figure('position', [100, 100, 1200, 1200]);
% figure/
hold on;

for t = 0:10
    % for t = [0:34:1256 1256]
    clf, hold on;
    % view(az, el);
    particles = dlmread([folder num2str(t)]);
    levelset = dlmread([folder 'ls' num2str(t)]);
    values = reshape(levelset, N, N)';

    [row, col] = find(particles(:, 5) == -1);
    negatives = particles(row, :);
    [row, col] = find(particles(:, 5) == +1);
    positives = particles(row, :);

    scatter(positives(:, 1), positives(:, 2), 'b', '.');
    scatter(negatives(:, 1), negatives(:, 2), 'r', '.');
    % quiver(particles(:, 1), particles(:, 2), particles(:, 2), -particles(:, 1));

    contour(X, Y, values, [0, 0], 'k-', 'linewidth', 2);
    title(t);
    pause;

    % pause;
    [az, el] = view();
end
