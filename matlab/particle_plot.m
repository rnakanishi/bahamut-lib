clear;

folder = '~/git/bahamut-lib/results/particles/2d/';
timesteps = 100;

az = 0;
el = 0;

figure('position', [100, 100, 1200, 600]);
% figure/
hold on;

for t = 0:10
    % for t = [0:34:1256 1256]
    clf, hold on;
    % view(az, el);
    particles = dlmread([folder num2str(t)]);

    scatter(particles(:, 1), particles(:, 2));
    quiver(particles(:, 1), particles(:, 2), particles(:, 2), -particles(:, 1));

    title(t);
    pause;

    % pause;
    [az, el] = view();
end
