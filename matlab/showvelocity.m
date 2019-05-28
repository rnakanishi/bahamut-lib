%N = 40;
x = linspace(0, 1, N);
[X, Y] = meshgrid(x, x);
hold on;
while true
    v = dlmread('../results/velocityField');
    clf;
    %quiver3(Y(:), X(:), Z(:), v(:,1), v(:,2), v(:,3), 5);
    quiver(Y(:), X(:), v(:,1), v(:,2), 5);
    axis equal tight;
    view(2);
    pause;
end