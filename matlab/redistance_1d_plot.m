clear;
N = 32;
timesteps = 300;

x = linspace(-5, 5, N);
X = x;
dx = x(2) - x(1);

az = 0;
el = 0;

figure('position', [100, 100, 1200, 600]);
% figure/
hold on;

for t = 0:timesteps
    clf, hold on;
    % view(az, el);
    levelset = dlmread(['~/git/bahamut-lib/results/redistance/1d/' num2str(t)]);

    values = levelset 

    set(gca, 'xtick', -5:dx:5);
    set(gca, 'ytick', -5:dx:5);
    % axis equal;
    axis([-5 5 -5 30]);
    grid on;
    title(num2str(t));
    plot(X, values);

    if (t == 0)
        pause
    else
        pause(0.2);
    end

    % pause;
end
