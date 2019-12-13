figure;
particles = dlmread("particles.txt");
% scatter3(particles(:, 1), particles(:, 2), particles(:, 3), '.');
% quiver3(particles(:, 1), particles(:, 2), particles(:, 3), particles(:, 4), particles(:, 5), particles(:, 6));
quiver(particles(:, 1), particles(:, 2), particles(:, 3), particles(:, 4));
% axis equal tight;

% grads = dlmread("gradients");
% quiver3(grads(:, 1), grads(:, 2), grads(:, 3), grads(:, 4), grads(:, 5), grads(:, 6), 5);
