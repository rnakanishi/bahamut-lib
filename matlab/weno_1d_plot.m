clear;
N = 25;
timesteps = 90;

x = linspace(0, 2 * pi, N);

figure, hold on;

for count = 0:timesteps - 1
    clf, hold on;
    weno = dlmread(['~/git/bahamut-lib/results/weno/1d/' num2str(count)]);

    analytic = weno(:, 1);
    weno = weno(:, 2);

    % axis([0 2*pi -2 2]);
    % set(anSurf, 'cData', zeros(1, N * N));
    plot(analytic)
    plot(weno)

    pause(0.1);
end
