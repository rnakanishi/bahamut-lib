addpath("~/matlab/wobj");

folders = '../results/particles/3d/enright/';
subfolders = ['semiL/'; 'weno/'; 'pls20/'; 'pls40/'; 'pls80/'];

for tech = 1:length(subfolders)

    folder = [folders deblank(subfolders(tech, :))];

    first = dlmread([folder 'ls' num2str(0)]);
    N = (length(first(:, 1)))^(1/3) + 1;

    for index = 0:181
        levelset = dlmread([folder 'ls' num2str(index)]);
        values = reshape(levelset(:, 4), N, N, N);
        values = permute(values, [2, 1, 3]);
        levelx = reshape(levelset(:, 2), N, N, N);
        levely = reshape(levelset(:, 1), N, N, N);
        levelz = reshape(levelset(:, 3), N, N, N);

        FV = isosurface(levelx, levely, levelz, values, 0);
        normals = isonormals(values, FV.vertices);

        L = sqrt(normals(:, 1).^2 + normals(:, 2).^2 + normals(:, 3).^2) + eps;
        normals(:, 1) = normals(:, 1) ./ L; normals(:, 2) = normals(:, 2) ./ L; normals(:, 3) = normals(:, 3) ./ L;
        face = FV.faces;
        vert = FV.vertices;

        y = zeros(length(FV.faces), 3);
        y(:, 1) = mean(reshape(normals(face(1:end, :), 1), 3, []));
        y(:, 2) = mean(reshape(normals(face(1:end, :), 2), 3, []));
        y(:, 3) = mean(reshape(normals(face(1:end, :), 3), 3, []));

        v1 = vert(face(:, 1), :) - vert(face(:, 2), :);
        v2 = vert(face(:, 1), :) - vert(face(:, 3), :);

        cNormal = cross(v1, v2);
        L = sqrt(cNormal(:, 1).^2 + cNormal(:, 2).^2 + cNormal(:, 3).^2) + eps;
        cNormal(:, 1) = cNormal(:, 1) ./ L; cNormal(:, 2) = cNormal(:, 2) ./ L; cNormal(:, 3) = cNormal(:, 3) ./ L;

        wrongNormals = sum(y .* cNormal, 2) < 0;
        face(wrongNormals, :) = face(wrongNormals, [1 3 2]);
        clear obj;
        obj.vertices = FV.vertices;
        % obj.vertices_normal = normals;
        obj.objects(1).type = 'f';
        obj.objects(1).data.vertices = face;
        % obj.objects(1).data.normal = face;
        objfilename = sprintf("%s%s%04d.obj", '../results/mesh/3d/enright/', ...
            deblank(subfolders(tech, :)), index);
        write_wobj(obj, objfilename);

        fprintf("Written %s\n", objfilename);
        % return
    end

end
