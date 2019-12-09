file = fopen("/home/rnakanishi/Documents/meshes/icosphere.sdf");

line = fgets(file);
dimension = sscanf(line, " %f %f %f ");

line = fgets(file);
origin = sscanf(line, " %f %f %f ");

line = fgets(file);
dx = sscanf(line, " %f ");

phi = fscanf(file, "%f");
phi = reshape(phi, dimension);

isosurface(phi, 0)
