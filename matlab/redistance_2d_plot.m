clear;
N = 20;
timesteps = 19;

x = linspace(-5, 5, N);
[X, Y] = meshgrid(x, x);
dx = x(2) - x(1);

az = 0;
el = 0;

figure('position', [100, 100, 1200, 600]);
% figure/
hold on;

for t = 0:timesteps
    clf, hold on;
    % view(az, el);
    levelset = dlmread(['~/git/bahamut-lib/results/redistance/2d/' num2str(t)]);

    values = reshape(levelset, N, N)';

    % anSurf = surf(X, Y, analytic);
    % subplot(121);
    % levelsurf = surf(X, Y, values);
    % axis([-5 5 -5 5 -5 5]);
    % axis equal;
    % title(num2str(t));

    % axis([-2*pi 2*pi -2*pi 2*pi]);
    % set(anSurf, 'cData', zeros(1, N * N));
    % set(weSurf, 'cData', ones(1, N * N));
    % xlabel('X');
    % ylabel('Y');

    % subplot(122);
    % hold on;
    contour(X, Y, values, -5:0.2:5);
    contour(X, Y, values, [0, 0], 'k', 'linewidth', 2);
    set(gca, 'xtick', -5:dx:5);
    set(gca, 'ytick', -5:dx:5);
    axis equal;
    axis([-5 5 -5 5]);
    grid on;
    title(num2str(t));

    if (t==0)
        pause
    else
        pause(0.01);
    end

    % pause;
    [az, el] = view();
end
