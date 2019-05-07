N = 40;
x = linspace(0, 1, N);

while true
    [X, Y, Z] = meshgrid(x, x, x);
    v = dlmread('divergent');
    v = reshape(v, N, N, N);
    v = permute(v, [1 3 2 ]);
    clf;
    zeroIds = find(v == 0);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    X(zeroIds) = [];
    Y(zeroIds) = [];
    Z(zeroIds) = [];
    v(zeroIds) = [];
    scatter3(X(:), Y(:), Z(:), 25, v(:), 'filled');
    axis([0 1 0 1 0 1]);
    colorbar;
    rotate3d on;
    pause;
end
