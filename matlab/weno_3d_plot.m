clear;
N = 32;
timesteps = 5;

x = linspace(-2 * pi, 2 * pi, N);
[X, Y] = meshgrid(x, x);

figure, hold on;
az = 0;
el = 90;

for count = 0:timesteps - 1
    clf, hold on;
    view(az, el);
    weno3d = dlmread(['~/git/bahamut-lib/results/weno/3d/' num2str(count)]);

    weno = reshape(weno3d, N, N, N);

    % anSurf = surf(X, Y, analytic);
    % weSurf = surf(X, Y, weno);
    isosurface(weno, 0.);
    axis equal;

    % axis([-2*pi 2*pi -2*pi 2*pi]);
    % set(anSurf, 'cData', zeros(1, N * N));
    % set(weSurf, 'cData', ones(1, N * N));
    xlabel('X');
    ylabel('Y');
    title(num2str(count))

    pause;
    % pause;
    [az, el] = view();
end
