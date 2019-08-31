clear;
N = 32;
timesteps = 300;

x = linspace(-2, 2, N);
[X, Y] = meshgrid(x, x);

az = 0;
el = 90;

figure(1, 'position', [100, 100, 1200, 600]);
hold on;

for count = 0:timesteps
    clf, hold on;
    view(az, el);
    levelset = dlmread(['~/git/bahamut-lib/results/redistance/2d/' num2str(count)]);

    values = reshape(levelset, N, N);

    % anSurf = surf(X, Y, analytic);
    % subplot(121);
    % levelsurf = surf(X, Y, values);
    % axis equal;

    % axis([-2*pi 2*pi -2*pi 2*pi]);
    % set(anSurf, 'cData', zeros(1, N * N));
    % set(weSurf, 'cData', ones(1, N * N));
    % xlabel('X');
    % ylabel('Y');
    % title(num2str(count));

    % subplot(122);
    contour(values, [0, 0]);
    axis equal;
    title(num2str(count));

    pause;
    % pause;
    [az, el] = view();
end
