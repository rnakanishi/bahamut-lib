clear;
N = 1000;
timesteps = 1000;

x = linspace(-5, 5, N);
X = x;
dx = x(2) - x(1);

az = 0;
el = 0;

figure('position', [100, 100, 1200, 600]);
% figure/
hold on;

for t = 150:timesteps
    clf, hold on;
    % view(az, el);
    levelset = dlmread(['~/git/bahamut-lib/results/cip/1d/' num2str(t)]);

    values = levelset;

    % set(gca, 'xtick', -5:dx:5);
    % set(gca, 'ytick', -5:dx:5);
    % axis equal;
    axis([-5 5 -2 7]);
    grid on;
    title(num2str(t));
    plot(X, values);

    print(['~/Documents/plots/cip/1d/' num2str(t) '.jpg'], '-djpg');
        pause(0.05)
    % pause;
end
