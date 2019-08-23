clear;
N = 50;
timesteps = 500;

x = linspace(-2 * pi, 2 * pi, N);
[X, Y] = meshgrid(x, x);

figure, hold on;
az = 0;
el = 90;

for count = 0:timesteps - 1
    clf, hold on;
    view(az, el);
    weno2d = dlmread(['~/git/bahamut-lib/results/weno/2d/' num2str(count)]);

    analytic = reshape(weno2d(:, 1), N, N);
    weno = reshape(weno2d(:, 2), N, N);

    % anSurf = surf(X, Y, analytic);
    % weSurf = surf(X, Y, weno);
    contour(weno, [0,0]);
    axis equal;

    % axis([-2*pi 2*pi -2*pi 2*pi]);
    % set(anSurf, 'cData', zeros(1, N * N));
    % set(weSurf, 'cData', ones(1, N * N));
    xlabel('X');
    ylabel('Y');
    title(num2str(count))

    pause(1/100);
    % pause;
    [az, el] = view();
end
