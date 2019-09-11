clear;

folder = '~/git/bahamut-lib/results/cip/2d/';
first = dlmread([folder num2str(0)]);
N = sqrt(length(first));
timesteps = 100;

x = linspace(-5, 5, N);
[X, Y] = meshgrid(x, x);
dx = x(2) - x(1);

az = 0;
el = 0;

figure('position', [100, 100, 1200, 600]);
% figure/
hold on;

for t = 0:5:1000
% for t = [0:34:1256 1256]
    clf, hold on;
    % view(az, el);
    levelset = dlmread([folder num2str(t)]);

    values = reshape(levelset, N, N)';

    % values = round(values * 1e10) / 1e10;

    % subplot(121);
    % levelsurf = surf(X, Y, values);
    % axis([-5 5 -5 5 -5 5]);
    % view(az, el);
    % axis equal;
    % title(num2str(t));

    % axis([-2*pi 2*pi -2*pi 2*pi]);
    % set(anSurf, 'cData', zeros(1, N * N));
    % set(weSurf, 'cData', ones(1, N * N));
    % xlabel('X');
    % ylabel('Y');

    % subplot(122);
    hold on;
    contour(X, Y, values, -1:0.2:1);
    contour(X, Y, values, [0, 0], 'k', 'linewidth', 2);
    set(gca, 'xtick', -5:dx:5);
    set(gca, 'ytick', -5:dx:5);
    axis equal;
    axis([-5 5 -5 5]);
    grid on;
    title(num2str(t));

    % if (t==0)
    pause
    % else
    % pause(0.05);
    % end

    % pause;
    [az, el] = view();
end
