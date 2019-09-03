N = 32;
seconds = 50;

x = linspace(-5, 5, N + 1);
dx = x(2) - x(1);

X = x;
X = (X(1:end - 1) + X(2:end)) / 2;
phi = (X - 0.4 * dx) .* (X + 6) / 2 + 1;
% phi = (sqrt(X.^2 + Y.^2) - 4);

phi0 = phi;

eps = dx;
dt = 0.9 * dx;
T = seconds / dt;

Interface = zeros(N, N);
interfaceL = [phi(1:end - 1) .* phi(2:end) 1];
interfaceR = [1 phi(2:end) .* phi(1:end - 1)];
Interface = (interfaceL <= 0) | (interfaceR <= 0) ;

S = phi0 ./ (phi0.^2 + eps^2);
Sbool = S;
Sbool(phi < 0) = -1;
Sbool(phi > 0) = 1;
Sbool(phi == 0) = 0;
S = Sbool;

centralDiffX = [0 (phi(3:end) - phi(1:end - 2)) 0];
D = dx * 2 * phi0 ./ abs(centralDiffX);

for t = 0:T
    clf;
    hold on;
    plot(X, phi);
    % axis equal;
    axis([-5 5 -5 30])
    set(gca, 'xtick', -5 + dx:dx:5 + dx);
    set(gca, 'ytick', -5 + dx:dx:5 + dx);
    grid on;

    % view(az, el);
    % surf(X, Y, phi);
    % axis([-5 5 -5 5 -5 5]);

    title([num2str(t) ' of ' num2str(T)])

    gradxdown = -fliplr(diff(fliplr(phi), 1, 2)) / dx;
    gradxup = diff(phi, 1, 2) / dx;

    a = [0 gradxdown];
    b = [gradxup 0];

    positive = max(abs(max(0 * a, a)), abs(min(0 * b, b))) - 1;
    negative = max(abs(min(0 * a, a)), abs(max(0 * b, b))) - 1;

    G = zeros(1, N);
    G(find(phi0 > 0)) = positive(find(phi0 > 0));
    G(find(phi0 < 0)) = negative(find(phi0 < 0));
    G(find(phi0 == 0)) = 0;

    phiNotInterface = phi - dt * S .* G;
    phiInterface = phi - dt / dx * (S .* abs(phi) - D);

    phi(Interface) = phiInterface(Interface);
    phi(~Interface) = phiNotInterface(~Interface);

    aa = max(0 * a, a);
    bb = min(0 * b, b);
    gradxp = max(aa, abs(bb));
    gradxp(gradxp == abs(bb)) = gradxp(gradxp == abs(bb)) .* sign(bb(gradxp == abs(bb)));

    aa = min(0 * a, a);
    bb = max(0 * b, b);
    gradxn = max(abs(a), b);
    gradxn(gradxn == abs(aa)) = gradxn(gradxn == abs(aa)) .* sign(aa(gradxn == abs(aa)));

    gradx = grady = zeros(1, N);
    gradx(phi0 > 0) = gradxp(phi0 > 0);
    gradx(phi0 < 0) = gradxn(phi0 < 0);

    quiver(X, 0*X, gradx, (0 * gradx)+0.2);

    if (t == 0)
        pause
    else
        pause(1/10);
    end

    [az, el] = view();
end
