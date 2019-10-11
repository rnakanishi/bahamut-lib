folder = "../results/";
files = [[folder "pls_log.csv"]; [folder "weno_log.csv"]; [folder "semiL_log.csv"]];
labels = ['Particle Levelset'; 'Weno advection'; 'Semi Lagrangean']

titles = ['Cell Advection'; 'Particle advection'];
titleIds = [1 6];

ndata = size(files, 1)
data = cell(1, ndata);

for nfile = 1:size(files, 1)
    _file = deblank(files(nfile, :));
    data{nfile} = csvread(_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Advection

fig = figure;
hold on;
title(titles(1, :))

plot(3 * data{1}(1, :));

for index = 2:ndata
    plot(2.25 * (1:length(data{index}(1, :))), data{index}(1, :));
end

legend(labels);
axis square;

imagename = sprintf("../results/images/graphs/cell_advection.jpg");
print(fig, imagename, "-S800,800");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Particle Advection x Number of particles
fig = figure;
hold on;
title(titles(2, :))

pAdvection = data{1}(6, :);
pCount = dlmread([folder 'particleCount_pls.txt']);

plotyy(1:length(pAdvection), pAdvection, 2.25 * (1:length(pCount)), pCount)
legend("Particle advection time", "Particle count");

axis square;

imagename = sprintf("../results/images/graphs/particle_advection.jpg");
print(fig, imagename, "-S800,800");
