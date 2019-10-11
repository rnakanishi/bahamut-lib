folder = "../results/";
files = [[folder "pls_log.csv"]; [folder "weno_log.csv"]; [folder "semiL_log.csv"]];
labels = ['Particle Levelset'; 'Weno advection']

titles = ['Cell Advection'; 'Particle advection'];
titleIds = [1 6];

ndata = size(files, 1)
data = cell(1, ndata);

for nfile = 1:size(files, 1)
    _file = deblank(files(nfile, :))
    data{nfile} = csvread(_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell Advection

figure;
hold on;
title(titles(1, :))

for index = 1:ndata
    plot(data{index}(1, :));
end

legend(labels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Particle Advection x Number of particles
figure;
hold on;
title(titles(2, :))

pAdvection = data{1}(6, :);
pCount = dlmread([folder 'particleCount_pls.txt']);

plotyy(1:length(pAdvection), pAdvection, 2.25 * (1:length(pCount)), pCount)
legend("Particle advection time", "Particle count");
