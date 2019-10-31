clear;

folder = '/home/rnakanishi/git/bahamut-lib/results/particles/3d/enright/cip/';
timesteps = 100;
first = dlmread([folder 'ls' num2str(0)]);
N = (length(first(:, 1)))^(1/3) + 1;

% x = linspace(0, 1, N);
x = linspace(-5, 5, N);
dx = x(2) - x(1);
[X, Y, Z] = meshgrid(x, x, x);

az = 0;
el = -90;

fig = figure('position', [10, 300, 1900, 1400]);
% figure/
hold on;

vis = [0 1 0 1 0 1];
% vis = [-5 5 -5 5 -5 5];

print("Starting loop")

for t = 2
    title([num2str(t) "..."]);
    drawnow;
    % for t = [100 150 200]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % particles = dlmread([folder 'p' num2str(t)]);
    % negatives = particles;

    % [row, col] = find(particles(:, 4) == -1);
    % negatives = particles(row, :);
    % [row, col] = find(particles(:, 4) == +1);
    % positives = particles(row, :);

    clf, hold on;
    % % subplot(122);
    % plot3(negatives(1:end, 1), negatives(1:end, 2), negatives(1:end, 3), 'markersize', 5, 'b.', 'LineStyle', 'none');
    % plot3(positives(1:end, 1), positives(1:end, 2), positives(1:end, 3), 'markersize', 0.5, 'r.', 'LineStyle', 'none' );
    % view(az, el);
    % axis(vis)
    % axis equal;
    % grid off;

    % quiver(particles(:, 1), particles(:, 2), particles(:, 2), -particles(:, 1));
    % quiver3(negatives(1:2:end, 1), negatives(1:2:end, 2), negatives(1:2:end, 3), 64 * negatives(1:2:end, 4), 64 * negatives(1:2:end, 5), 64 * negatives(1:2:end, 6));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    levelset = dlmread([folder 'ls' num2str(t)]);
    values = reshape(levelset(:, 4), N, N, N);
    values = permute(values, [2, 1, 3]);
    levelx = reshape(levelset(:, 2), N, N, N);
    levely = reshape(levelset(:, 1), N, N, N);
    levelz = reshape(levelset(:, 3), N, N, N);

    % % subplot(121);
    clf;
    hold on;
    view(az, el);
    [face, vert] = isosurface(levelx, levely, levelz, values, 0);
    % p = patch("Faces", face, "Vertices", vert, "EdgeColor", "none", "facealpha", 0.7);
    p = trisurf(face, vert(:, 1), vert(:, 2), vert(:, 3));
    pbaspect([1 1 1]);
    isonormals(levelx, levely, levelz, values, p);
    set(p, "facecolor", "cyan", "edgecolor", "none", "FaceLighting", "gouraud");
    material dull
    light("Position", [3 -.5 .5], "color", [0.95 0.95 0.95]);
    light("Position", [-0.25 1 -.5], "color", [0.6 0.6 0.6]);

    % % light("Position", [0 0 10], "color", [0.95 0.95 0.95]);
    % % light("Position", [-5 -5 0], "color", [0.6 0.6 0.6]);

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
    axis off;
    rotate3d on;
    pause(0.1);
    % imagename = sprintf("../results/images/3d/pls_deform_lagrangean/%04d.jpg", t);
    % print(fig, imagename, "-S1150,1000");
    % fprintf("Written: %s\n", imagename);

    % pause;
    [az, el] = view();
end

% levelset = dlmread([folder 'ls' num2str(0)]);
% values = reshape(levelset(:, 3), N, N)';
% levelx = reshape(levelset(:, 1), N, N)';
% levely = reshape(levelset(:, 2), N, N)';
% hold on;
% contour(levelx, levely, values, [0, 0], 'k:', 'linewidth', 3);
