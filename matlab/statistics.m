addpath("~/matlab/cbrewer/")

folder = "../results/statistics/enright128/";
files = [
% folder "pls_log20.csv";
    folder "pls64_log.csv";
% folder "pls_log80.csv";
    folder "weno_log.csv";
    folder "semiL_log.csv";
    folder "cip_log.csv"
    ];
labels = [
% 'Particle Levelset 20';
    'Particle Levelset 64';
% 'Particle Levelset 80';
    'Weno advection';
    'Semi Lagrangean';
    'CIP'
    ];

titles = ['Cell Advection'; 'Cell advection x Surface cells'; 'Particle advection'; 'Total'];

nfiles = size(files, 1);
data = cell(1, nfiles);

colors = cbrewer('qual', 'Set1', 5);

for nfile = 1:nfiles
    _file = deblank(files(nfile, :));
    data{nfile} = csvread(_file);
end

fig = figure('position', [10, 300, 1900, 1400]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total average time

subplot(221); % figure;
hold on;
title(deblank(titles(4, :)), 'fontsize', 14)

for index = 1:nfiles
    totaltime = data{index}(2:end, :);
    totaltime = sum(totaltime, 2) - totaltime(:, end - 1);
    h = plot((1:length(totaltime)), totaltime);
    set(h, 'color', colors(index, :));
    set(h, 'linewidth', 1.5);
end

xlabel('Timesteps', 'fontsize', 14)
ylabel('Seconds', 'fontsize', 14)
h = legend(labels, 'location', 'northeast');
set(h, 'fontsize', 12)

axis square;

% imagename = sprintf("../results/images/graphs/total.jpg");
% print(fig, imagename, "-S800,800");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Advection

subplot(222); % figure;
hold on;
title(deblank(titles(1, :)), 'fontsize', 14)

for index = 1:nfiles
    h = plot((1:length(data{index}(2:end, 1))), data{index}(2:end, 1));
    set(h, 'color', colors(index, :));
    set(h, 'linewidth', 1.5);
end

xlabel('Timesteps', 'fontsize', 14)
ylabel('Seconds', 'fontsize', 14)
h = legend(labels, 'location', 'northeast');
set(h, 'fontsize', 12)

axis square;

% imagename = sprintf("../results/images/graphs/cell_advection.jpg");
% print(fig, imagename, "-S800,800");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell  Advection x Number of cells
subplot(223); % figure;
hold on;
title(deblank(titles(2, :)), 'fontsize', 14)

cells = ['weno'; 'semiL'; 'pls'];

for index = 1:nfiles

    cAdvection = data{index}(2:end, 1);
    cCount = dlmread([folder  deblank(cells(index, :)) 'cellCount.txt']);

    [ax, h1, h2] = plotyy(1:length(cAdvection), cAdvection, (1:length(cCount)), cCount);
    set(h1, 'color', colors(index, :));
    set(h2, 'color', colors(index, :));
    set(h2, 'linestyle', '--');
    set([h1, h2], 'linewidth', 1.5);
end

xlabel('Timesteps', 'fontsize', 14)
ylabel(ax(1), 'Seconds', 'fontsize', 14)
ylabel(ax(2), 'Number of surface cells', 'fontsize', 14)
set(ax, 'ycolor', 'k')
h = legend([h1, h2], 'Advection', 'NCells')
set(h, 'fontsize', 12)

% legend("Advection 20", "Advection 40", "Advection 80",
%    "Particle count 20", "Particle count 40", "Particle count 80",
%    'location', 'northeast');

axis square;

% imagename = sprintf("../results/images/graphs/particle_advection.jpg");
% print(fig, imagename, "-S800,800");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Particle Advection x Number of particles
subplot(224); % figure;
hold on;
title(deblank(titles(3, :)), 'fontsize', 14)

particles = ['20'; '40'];

for index = 1:3

    pAdvection = data{index}(2:end, 6);
    pCount = dlmread([folder 'particleCount' particles(index, :) '.txt']);

    [~, ids] = unique(pCount);
    pCount = pCount(sort(ids));

    [ax, h1, h2] = plotyy(1:length(pAdvection), pAdvection, (1:length(pCount)), pCount);
    set([h1, h2], 'color', colors(index, :));
    set(h2, 'linestyle', '--');
    set([h1, h2], 'linewidth', 1.5);
end

xlabel('Timesteps', 'fontsize', 14)
ylabel(ax(1), 'Seconds', 'fontsize', 14)
ylabel(ax(2), 'Number of particles', 'fontsize', 14)
set(ax, 'ycolor', 'k')
h = legend("Advection 20", "Advection 40", "Advection 80", ...
    "Particle count 20", "Particle count 40", "Particle count 80", ...
    'location', 'northeast');
set(h, 'fontsize', 12)

axis square;

% imagename = sprintf("../results/images/graphs/particle_advection.jpg");
% print(fig, imagename, "-S800,800");

imagename = sprintf("../results/images/graphs/statistics.jpg");
print(fig, imagename, "-S2000,2000");
