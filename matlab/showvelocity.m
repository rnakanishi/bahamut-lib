N = 40;
x = linspace(0, 1, N);
[X, Y, Z] = meshgrid(x, x, x);

while true
    v = dlmread('velocityField');
    clf;
    quiver3(Y(:), X(:), Z(:), v(:,1), v(:,2), v(:,3), 5);
    axis equal tight;
    view(2);
    pause;
end