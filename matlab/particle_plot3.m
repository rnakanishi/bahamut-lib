clear;

folder = '~/git/bahamut-lib/results/particles/3d/pls_deform/';
timesteps = 100;
first = dlmread([folder 'ls' num2str(0)]);
N = (length(first(:, 1)))^(1/3) + 1;

x = linspace(0, 1, N);
% x = linspace(-5, 5, N);
dx = x(2) - x(1);
[X, Y, Z] = meshgrid(x, x, x);

az = 0;
el = 0;

fig = figure('position', [350, 150, 1500, 1500]);
% figure/
hold on;

vis = [0 1 0 1 0 1];
% vis = [-5 5 -5 5 -5 5];

for t = 40:1:90
    % for t = [100 150 200]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    particles = dlmread([folder num2str(t)]);

    [row, col] = find(particles(:, 7) == -1);
    negatives = particles(row, :);
    [row, col] = find(particles(:, 7) == +1);
    positives = particles(row, :);

    clf, hold on;

    plot3(negatives(1:2:end, 1), negatives(1:2:end, 2), negatives(1:2:end, 3), 'markersize', 0.5, 'b.', 'LineStyle', 'none');
    % plot3(positives(1:10:end, 1), positives(1:10:end, 2), positives(1:10:end, 3), 'markersize', 0.5, 'r.', 'LineStyle', 'none' );

    % quiver(particles(:, 1), particles(:, 2), particles(:, 2), -particles(:, 1));
    % quiver3(negatives(1:2:end, 1), negatives(1:2:end, 2), negatives(1:2:end, 3), 64 * negatives(1:2:end, 4), 64 * negatives(1:2:end, 5), 64 * negatives(1:2:end, 6));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    levelset = dlmread([folder 'ls' num2str(t)]);
    values = reshape(levelset(:, 4), N, N, N);
    values = permute(values, [2,1,3]);
    levelx = reshape(levelset(:, 2), N, N, N);
    levely = reshape(levelset(:, 1), N, N, N);
    levelz = reshape(levelset(:, 3), N, N, N);
    
    clf; 
    hold on;
    view(az, el);
    [face, vert] = isosurface(levelx, levely, levelz, values, 0);
    p = patch("Faces", face, "Vertices", vert, "EdgeColor", "none");
    pbaspect([1 1 1]);
    isonormals(levelx, levely, levelz, values, p);
    set(p, "FaceColor", "cyan", "FaceLighting", "gouraud");
    light("Position", [1 1 5]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % grad = dlmread([folder 'v' num2str(t)]);
    % scale = 5;
    % gradx = reshape(grad(:, 4), N, N, N);
    % grady = reshape(grad(:, 5), N, N, N);
    % gradz = reshape(grad(:, 6), N, N, N);
    % positionx = reshape(grad(:, 1), N, N, N);
    % positiony = reshape(grad(:, 2), N, N, N);
    % positionz = reshape(grad(:, 3), N, N, N);
    % hold on;
    % quiver3(positionx, positiony, positionz, scale * gradx, scale * grady, scale * gradz, 'k');

    % [V,U , W] = gradient(values);
    % quiver3(positionx, positiony, positionz, scale * U, scale * V, scale * W, 'r');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % grad = dlmread([folder 'v' num2str(t)]);
    % scale = 2.5;
    % gradx = reshape(grad(:, 3), N, N)';
    % grady = reshape(grad(:, 4), N, N)';
    % positionx = reshape(grad(:, 1), N, N)';
    % positiony = reshape(grad(:, 2), N, N)';
    % hold on;
    % quiver(positionx, positiony, scale * gradx, scale * grady, 'y');

    % [U, V] = gradient(values);
    % quiver(positionx, positiony, scale*U, scale*V, 'r');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axis equal
    title(t);
    axis(vis);
    pause(0);
    axis off;
    % print(fig, ["../results/images/3d/pls_rigid/" int2str(t) ".jpg"], "-S1750,1800")

    % pause;
    [az, el] = view();
end

levelset = dlmread([folder 'ls' num2str(0)]);
values = reshape(levelset(:, 3), N, N)';
levelx = reshape(levelset(:, 1), N, N)';
levely = reshape(levelset(:, 2), N, N)';
hold on;
contour(levelx, levely, values, [0, 0], 'k:', 'linewidth', 3);
