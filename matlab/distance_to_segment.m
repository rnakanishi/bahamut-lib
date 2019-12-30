function [d t] = distance_to_segment(origin, ending, target)
    figure, hold on, axis equal 
    points = [origin; ending; target];
    scatter(points(:, 1), points(:, 2))

    vec = ending - origin;
    tvec = target - origin;
    t = tvec * vec' / (norm(vec)^2);

    quiver(origin(1), origin(2), vec(1), vec(2), 'b')
    quiver(origin(1), origin(2), tvec(1), tvec(2), 'r')

    t = min(1, max(0, t));
    proj = origin + t * vec;
    dvec = proj - target;

    quiver(target(1), target(2), dvec(1), dvec(2), 'g')

    display(norm(origin - target))

    scatter(proj(1), proj(2), '*')

    tvec = tvec / norm(tvec);
    vec = vec / norm(vec);

    d = norm(dvec);

    if (tvec(1) * vec(2) - tvec(2) * vec(1) < 0)
        d = -d;
    end
