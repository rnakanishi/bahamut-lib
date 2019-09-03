N = 20;
seconds = 10;

x = linspace(-5, 5, N + 1);
[X, Y] = meshgrid(x, x);

X = (X(:, 1:end - 1) + X(:, 2:end)) / 2;
Y = (Y(:, 1:end - 1) + Y(:, 2:end)) / 2;
X = (X(1:end - 1, :) + X(2:end, :)) / 2;
Y = (Y(1:end - 1, :) + Y(2:end, :)) / 2;

f = 0.5 + (X - 3.5).^2 + (Y - 2).^2;
A = 4;
B = 2;
phi = f .* (sqrt(X.^2 ./ A^2 + Y.^2 ./ B^2) - 1);
% phi = (sqrt(X.^2 + Y.^2) - 4);

phi0 = phi;

dx = x(2) - x(1);
eps = dx;
dt = 0.5 * dx;
T = seconds / dt;

Interface = zeros(N, N);
interfaceL = [phi(:, 1:end - 1) .* phi(:, 2:end) ones(N, 1)];
interfaceR = [ones(N, 1) phi(:, 2:end) .* phi(:, 1:end - 1)];
interfaceT = [phi(1:end - 1, :) .* phi(2:end, :); ones(1, N)];
interfaceB = [ones(1, N); phi(2:end, :) .* phi(1:end - 1, :)];
Interface = (interfaceL <= 0) | (interfaceR <= 0) | (interfaceT <= 0) | (interfaceB <= 0);

S = phi0 ./ (phi0.^2 + eps^2);
Sbool = S;
Sbool(phi < 0) = -1;
Sbool(phi > 0) = 1;
Sbool(phi == 0) = 0;
S = Sbool;

centralDiffX = [zeros(N, 1) (phi(:, 3:end) - phi(:, 1:end - 2)) zeros(N, 1)];
centralDiffY = [zeros(1, N); (phi(3:end, :) - phi(1:end - 2, :)); zeros(1, N)];
D = dx * 2 * phi0 ./ sqrt(centralDiffX.^2 + centralDiffY.^2);

% S(Interface) = D(Interface)/ dx;

az = 0;
el = 0;

for t = 0:T + 1
    clf;
    hold on;
    contour(X, Y, phi, -5:0.2:5);
    contour(X, Y, phi, [0, 0], 'k', 'linewidth', 2);
    % axis equal;
    axis([-5 5 -5 5])
    set(gca, 'xtick', -5 + dx:dx:5 + dx);
    set(gca, 'ytick', -5 + dx:dx:5 + dx);
    grid on;

    % view(az, el);
    % surf(X, Y, phi);
    % axis([-5 5 -5 5 -5 5]);

    title([num2str(t) ' of ' num2str(T)])

    gradxdown = -fliplr(diff(fliplr(phi), 1, 2)) / dx;
    gradydown = -flipud(diff(flipud(phi), 1, 1)) / dx;
    gradxup = diff(phi, 1, 2) / dx;
    gradyup = diff(phi, 1, 1) / dx;

    a = [zeros(N, 1) gradxdown];
    b = [gradxup zeros(N, 1)];
    c = [zeros(1, N); gradydown];
    d = [gradyup; zeros(1, N)];

    positive = sqrt(max(cat(3, max(0 * a, a).^2, min(0 * b, b).^2), [], 3) + ...
        max(cat(3, max(0 * c, c).^2, min(0 * d, d).^2), [], 3)) - 1;
    negative = sqrt(max(cat(3, min(0 * a, a).^2, max(0 * b, b).^2), [], 3) + ...
        max(cat(3, min(0 * c, c).^2, max(0 * d, d).^2), [], 3)) - 1;

    G = zeros(N, N);
    G(find(phi0 > 0)) = positive(find(phi0 > 0));
    G(find(phi0 < 0)) = negative(find(phi0 < 0));
    G(find(phi0 == 0)) = 0;

    G = reshape(G, N, N);

    phiNotInterface = phi - dt * S .* G;
    phiInterface = phi - dt / dx * (S .* abs(phi) - D);

    phi(Interface) = phiInterface(Interface);
    phi(~Interface) = phiNotInterface(~Interface);

    aa = max(0 * a, a);
    bb = min(0 * b, b);
    cc = max(0 * c, c);
    dd = min(0 * d, d);
    gradxp = max(aa, abs(bb));
    gradyp = max(cc, abs(dd));
    gradxp(gradxp == abs(bb)) = gradxp(gradxp == abs(bb)) .* sign(bb(gradxp == abs(bb)));
    gradyp(gradyp == abs(dd)) = gradyp(gradyp == abs(dd)) .* sign(dd(gradyp == abs(dd)));

    aa = min(0 * a, a);
    bb = max(0 * b, b);
    cc = min(0 * c, c);
    dd = max(0 * d, d);
    gradxn = max(abs(a), b);
    gradyn = max(abs(c), d);
    gradxn(gradxn == abs(aa)) = gradxn(gradxn == abs(aa)) .* sign(aa(gradxn == abs(aa)));
    gradyn(gradyn == abs(cc)) = gradyn(gradyn == abs(cc)) .* sign(cc(gradyn == abs(cc)));

    gradx = grady = zeros(N, N);
    gradx(phi0 > 0) = gradxp(phi0 > 0);
    gradx(phi0 < 0) = gradxn(phi0 < 0);
    grady(phi0 > 0) = gradyp(phi0 > 0);
    grady(phi0 < 0) = gradyn(phi0 < 0);

    quiver(X, Y, gradx, grady);

    % if (t == 0)
    pause
    % else
    %     pause(1/20);
    % end

    [az, el] = view();
end
